import streamlit as st
import re
import pandas as pd
import io
import requests
import difflib
import xml.etree.ElementTree as ET
from stmol import showmol
import py3Dmol

# ==========================================
# 1. 网页全局设置
# ==========================================
st.set_page_config(page_title="工业级抗体生信大屏 V13.1", page_icon="🧬", layout="wide")

st.title("🧬 工业级抗体生信、CMC 与 IP 联合防御系统 (V13.1 完整版)")
st.info("💡 全栈闭环工作流：集成 HMM/Regex 双核提取、高阶成药性评估、多维聚类、Excel 序列清洗、WIPO 专利破译 (ST.26/ST.25)，以及本地 FTO 专利侵权排查雷达。")

# ==========================================
# 2. 深度 CMC 评估引擎
# ==========================================
def calculate_pi(seq):
    seq = seq.upper()
    pKa_acidic = {'D': 3.9, 'E': 4.1, 'C': 8.5, 'Y': 10.1}
    pKa_basic = {'K': 10.8, 'R': 12.5, 'H': 6.5}
    def net_charge(pH):
        charge = 1.0 / (1.0 + 10**(pH - 8.0)) - 1.0 / (1.0 + 10**(3.1 - pH))
        for aa, pka in pKa_acidic.items(): charge -= seq.count(aa) / (1.0 + 10**(pka - pH))
        for aa, pka in pKa_basic.items(): charge += seq.count(aa) / (1.0 + 10**(pH - pka))
        return charge
    pH_min, pH_max = 0.0, 14.0
    for _ in range(100):
        pH_mid = (pH_min + pH_max) / 2
        if net_charge(pH_mid) > 0: pH_min = pH_mid
        else: pH_max = pH_mid
    return round(pH_mid, 2)

def calculate_gravy(seq):
    if not seq or seq == "未识别": return 0.0
    kd_scale = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
    return round(sum(kd_scale.get(aa, 0.0) for aa in seq.upper()) / len(seq), 3)

def detect_unpaired_cysteine(seq):
    cys_positions = [i+1 for i, aa in enumerate(seq.upper()) if aa == 'C']
    count = len(cys_positions)
    if count % 2 != 0: return f"🚨 高危: 奇数({count})个 Cys @{cys_positions}"
    elif count > 2 and count % 2 == 0: return f"⚠️ 警告: 额外配对({count})个 Cys @{cys_positions}"
    return "✅ 正常 (2个 Cys)"

def guess_germline(seq):
    seq = seq.upper()[:50]
    if re.search(r'GGGSVQ', seq) or re.search(r'W[FY]RQAPGKERE', seq): return "Camelid VHH (纳米抗体)"
    if re.search(r'SGGGLVQ', seq): return "Human IGHV3 (高人源化)"
    if re.search(r'SGAEVKKPG', seq): return "Human IGHV1/5 (高人源化)"
    if re.search(r'SGSELKKPG', seq): return "Human IGHV7 (高人源化)"
    if re.search(r'SGPGLVKPSG', seq): return "Human IGHV4 (高人源化)"
    if re.search(r'SGPEVKKPG', seq): return "Human IGHV2 (高人源化)"
    if re.search(r'QQSG[AP]E[LV]V', seq) or re.search(r'QQSDA', seq) or re.search(r'GSLKLS', seq): return "Murine IGHV (鼠源)"
    if re.search(r'VQL[VQE]QSG', seq) or re.search(r'VQL[LVE]ESG', seq): return "IGHV (亚族未定)"
    if re.search(r'SP[SS][SF]LSASVG', seq): return "Human IGKV1/3 (高人源化)"
    if re.search(r'SPLSLPVTPG', seq): return "Human IGKV2 (高人源化)"
    if re.search(r'SP[DS]SLA[VS]SLG', seq): return "Human IGKV4 (高人源化)"
    if re.search(r'SP[SA]YLAASP', seq) or re.search(r'FMSTSVG', seq): return "Murine IGKV (鼠源)"
    if re.search(r'[DE][IV][VQAM][ML]TQS', seq): return "IGKV (亚族未定)"
    if re.search(r'QPPS[AS]SG', seq) or re.search(r'Q[PP]SVS[VAS]P', seq): return "Human IGLV (Lambda)"
    if re.search(r'LTQP', seq): return "IGLV (Lambda未定)"
    if seq.startswith('Q') or seq.startswith('E'): return "疑似重链 (高度变异)"
    if seq.startswith('D') or seq.startswith('A'): return "疑似轻链 (高度变异)"
    return "未知架构"

def get_region_finder(seq, cdrs, domain_type):
    if "Fc" in domain_type: return lambda i: "Fc区"
    c1, c2, c3 = cdrs.get("CDR1",""), cdrs.get("CDR2",""), cdrs.get("CDR3","")
    idx1 = seq.find(c1) if c1 != "未识别" else -1
    idx2 = seq.find(c2, max(0, idx1)) if c2 != "未识别" else -1
    idx3 = seq.find(c3, max(0, idx2)) if c3 != "未识别" else -1
    
    def region_of(i):
        if idx1 != -1 and i < idx1: return "FR1"
        if idx1 != -1 and i < idx1 + len(c1): return "CDR1"
        if idx1 != -1 and idx2 != -1 and i < idx2: return "FR2"
        if idx2 != -1 and i < idx2 + len(c2): return "CDR2"
        if idx2 != -1 and idx3 != -1 and i < idx3: return "FR3"
        if idx3 != -1 and i < idx3 + len(c3): return "CDR3"
        if idx3 != -1 and i >= idx3 + len(c3): return "FR4"
        return "可变区"
    return region_of

def detect_ptms_detailed(seq, cdrs, domain_type):
    region_finder = get_region_finder(seq, cdrs, domain_type)
    ptm_rules = {"N-糖基化": r"N[^P][ST]", "脱氨基": r"N[GSN]", "异构化": r"D[GS]", "酸断裂": r"DP", "氧化": r"M"}
    found_ptms = []
    for ptm_name, pattern in ptm_rules.items():
        for match in re.finditer(pattern, seq):
            found_ptms.append(f"[{region_finder(match.start())}] {ptm_name}({match.group()}) @{match.start()+1}")
    found_ptms.sort(key=lambda x: int(re.search(r'@(\d+)', x).group(1)) if re.search(r'@(\d+)', x) else 0)
    return " | ".join(found_ptms) if found_ptms else "✅ 无常见高危 PTM"

# ==========================================
# 3. 提取引擎 (API 云端引擎 vs 本地正则)
# ==========================================
def extract_cdrs_via_api(seq):
    """云端 API 请求器 (带自动降级 Fallback)"""
    api_url = "https://api.antibody-informatics.org/v1/anarci/annotate"
    payload = {"sequence": seq, "scheme": "imgt"}
    try:
        response = requests.post(api_url, json=payload, headers={"Content-Type": "application/json"}, timeout=3)
        if response.status_code == 200:
            data = response.json()
            return { "CDR1": data.get("CDR1", "未识别"), "CDR2": data.get("CDR2", "未识别"), "CDR3": data.get("CDR3", "未识别") }
        else: raise Exception("API 异常")
    except Exception:
        if "M" in seq or "L" in seq:
             return extract_vh_cdrs_regex(seq) if seq.startswith("E") or seq.startswith("Q") else extract_vl_cdrs_regex(seq)
        return extract_vh_cdrs_regex(seq)

def extract_vh_cdrs_regex(vh_seq):
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别"}
    cdr3_match = re.search(r"Y[YFCA]C[A-Z]{1,3}(.*?)W[GS][A-Z][GTSVI]", vh_seq)
    if cdr3_match: cdrs["CDR3"] = cdr3_match.group(1)
    cdr1_match = re.search(r"C[A-Z]{2,6}(.{5,16})W[A-Z][RQK]", vh_seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(1)
    cdr2_match = re.search(r"(?:EW[IVMASTL][A-Z]|KWM[A-Z]|REG[VLIA][A-Z]|RWV[A-Z])(.{8,30}?)[RKQ][VFSTIAM][TVIAMFSC][A-Z]?", vh_seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    return cdrs

def extract_vl_cdrs_regex(vl_seq):
    cdrs = {"CDR1": "未识别", "CDR2": "未识别", "CDR3": "未识别"}
    cdr3_match = re.search(r"Y[YFCA]C(.*?)(?:F[GSA][A-Z][GTV]|FGC)", vl_seq)
    if cdr3_match: cdrs["CDR3"] = cdr3_match.group(1)
    cdr1_match = re.search(r"C(.{8,18})W[YFL]", vl_seq)
    if cdr1_match: cdrs["CDR1"] = cdr1_match.group(1)
    cdr2_match = re.search(r"[ILVM][A-Z]([A-Z]{7})G[A-Z]P", vl_seq)
    if not cdr2_match: cdr2_match = re.search(r"W[YFL].{10,22}?([A-Z]{7})G[A-Z]{1,2}[RFS]", vl_seq)
    if cdr2_match: cdrs["CDR2"] = cdr2_match.group(1)
    return cdrs

def parse_fasta(text):
    sequences = {}
    if ">" not in text:
        sequences["未命名序列_1"] = re.sub(r'\s+', '', text).upper()
        return sequences
    for part in text.split(">"):
        if not part.strip(): continue
        lines = part.strip().split("\n")
        name, seq = lines[0].strip(), "".join(lines[1:]).replace(" ", "").upper()
        if name and seq: sequences[name] = seq
    return sequences
# ==========================================
# 3.5 3D 结构建模引擎 (ESMFold API)
# ==========================================
def render_3d_structure(pdb_string):
    """在线渲染 3D 结构"""
    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb_string, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    showmol(view, height=500, width=800)

@st.cache_data
def fetch_esm_fold_pdb(sequence):
    """调用 ESMFold 获取原子坐标"""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        response = requests.post(url, data=sequence, timeout=15)
        if response.status_code == 200: return response.text
    except: pass
    return None
# ==========================================
# 4. 主线战区：核心抗体解析与 CMC 聚类
# ==========================================
st.markdown("### 🔬 核心引擎：序列多维解析与 CMC 智能聚类")
engine_choice = st.radio("⚙️ 请选择底层提取引擎：", 
                         ("🚀 本地安全局域网版 (Regex 引擎) - 零延迟，极度保护 IP，适合内网跑批", 
                          "☁️ 外部云端对齐版 (RESTful API) - 智能发送至云端服务器打分 (带超时自动降级)"))

raw_input = st.text_area("📥 请在此粘贴需要评估的候选序列 (支持 FASTA 格式):", height=200, key="main_input")

if st.button("🚀 启动深度解析与聚类计算", type="primary"):
    use_api = "API" in engine_choice
    
    if raw_input:
        seq_dict = parse_fasta(raw_input)
        st.success(f"✅ 检测到 {len(seq_dict)} 条输入序列，正在使用 {'云端 API 引擎' if use_api else '本地高并发 Regex'} 进行处理...")
        
        all_results = []
        vh_pattern = re.compile(r"([EQ].{100,135}VTVSS)")
        vl_pattern = re.compile(r"([DEQA].{95,125}(?:VEIK|LEIK|TVLG|VTVL|FGC))")
        fc_pattern = re.compile(r"(CPPCP.*?LSPGK)")
        
        progress_bar = st.progress(0)
        total_seqs = len(seq_dict)
        
        for idx, (seq_name, clean_seq) in enumerate(seq_dict.items()):
            for i, vh in enumerate(vh_pattern.findall(clean_seq)):
                cdrs = extract_cdrs_via_api(vh) if use_api else extract_vh_cdrs_regex(vh)
                comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
                all_results.append({
                    "序列名称 (FASTA ID)": seq_name, "结构域": f"VH_{i+1}", "抽象结构域": "重链/纳米抗体",
                    "推测 Germline": guess_germline(vh), "孤立Cys 雷达": detect_unpaired_cysteine(vh),
                    "CDR_GRAVY": calculate_gravy(comb_cdr), "pI": calculate_pi(vh),
                    "PTM 空间定位分析": detect_ptms_detailed(vh, cdrs, "VH"),
                    "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": vh
                })
            for i, vl in enumerate(vl_pattern.findall(clean_seq)):
                cdrs = extract_cdrs_via_api(vl) if use_api else extract_vl_cdrs_regex(vl)
                comb_cdr = (cdrs["CDR1"] + cdrs["CDR2"] + cdrs["CDR3"]).replace("未识别", "")
                all_results.append({
                    "序列名称 (FASTA ID)": seq_name, "结构域": f"VL_{i+1}", "抽象结构域": "轻链",
                    "推测 Germline": guess_germline(vl), "孤立Cys 雷达": detect_unpaired_cysteine(vl),
                    "CDR_GRAVY": calculate_gravy(comb_cdr), "pI": calculate_pi(vl),
                    "PTM 空间定位分析": detect_ptms_detailed(vl, cdrs, "VL"),
                    "CDR1": cdrs["CDR1"], "CDR2": cdrs["CDR2"], "CDR3": cdrs["CDR3"], "完整序列": vl
                })
            for i, fc in enumerate(fc_pattern.findall(clean_seq)):
                all_results.append({
                    "序列名称 (FASTA ID)": seq_name, "结构域": f"Fc_{i+1}", "抽象结构域": "Fc区",
                    "推测 Germline": "IgG Fc", "孤立Cys 雷达": "-", "CDR_GRAVY": "-",
                    "pI": calculate_pi(fc), "PTM 空间定位分析": detect_ptms_detailed(fc, {}, "Fc"),
                    "CDR1": "-", "CDR2": "-", "CDR3": "-", "完整序列": fc
                })
            progress_bar.progress((idx + 1) / total_seqs)
        
        if all_results:
            df = pd.DataFrame(all_results)
            df_v = df[df['抽象结构域'].isin(['重链/纳米抗体', '轻链'])]
            cluster_v = df_v.groupby('完整序列').agg(
                抽象结构域=('抽象结构域', 'first'), 包含相同数量=('序列名称 (FASTA ID)', 'count'),
                来源分子名单=('序列名称 (FASTA ID)', lambda x: ', '.join(x.unique())),
                孤立Cys雷达=('孤立Cys 雷达', 'first'), CDR_GRAVY=('CDR_GRAVY', 'first'),
                CDR1=('CDR1', 'first'), CDR2=('CDR2', 'first'), CDR3=('CDR3', 'first'), Germline=('推测 Germline', 'first')
            ).reset_index().sort_values(by=['抽象结构域', '包含相同数量'], ascending=[True, False])
            
            df_valid_cdr3 = df_v[~df_v['CDR3'].str.contains('未识别|失败', na=False)]
            cluster_cdr3 = df_valid_cdr3.groupby('CDR3').agg(
                抽象结构域=('抽象结构域', 'first'), 共享该CDR3数量=('序列名称 (FASTA ID)', 'count'),
                来源分子名单=('序列名称 (FASTA ID)', lambda x: ', '.join(x.unique())), 代表完整序列=('完整序列', 'first')
            ).reset_index().sort_values(by=['抽象结构域', '共享该CDR3数量'], ascending=[True, False])
            
            st.markdown("#### 📊 CMC 解析与聚类总表")
            st.dataframe(df, use_container_width=True)
            
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                df.drop(columns=['抽象结构域']).to_excel(writer, index=False, sheet_name='Total_Report_总表')
                if not cluster_v.empty: cluster_v.to_excel(writer, index=False, sheet_name='Unique_V_Regions')
                if not cluster_cdr3.empty: cluster_cdr3.to_excel(writer, index=False, sheet_name='Unique_CDR3')
            
            st.download_button("📥 导出综合成药性与聚类报告 (.xlsx)", data=buffer.getvalue(), file_name="Antibody_CMC_Report.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="primary")
        else: st.warning("未能识别出标准片段。")
    else: st.error("请先输入序列！")
        # ... (上面是导出 Excel 的代码) ...
        else: st.warning("未能识别出标准片段。")
    else: st.error("请先输入序列！")

# 👇 从这里开始插入 3D 实验室前端代码 👇
# ==========================================
# 4.5 🧊 3D 结构实验室：即时建模与可视化
# ==========================================
st.markdown("---")
st.markdown("### 🧊 3D 结构实验室：即时建模与可视化")

# 只有当主战区跑完，生成了 all_results 且不为空时，才显示 3D 模块
if 'all_results' in locals() and all_results:
    st.info("💡 选中上方解析出的任意一条重链或轻链，一键调用 ESMFold 预测三维原子结构。")
    
    # 提取所有成功解析的序列名字作为下拉菜单选项
    mol_options = [f"{r['序列名称 (FASTA ID)']} ({r['结构域']})" for r in all_results if r['完整序列'] != "-" and len(r['完整序列']) > 20]
    
    if mol_options:
        selected_mol_name = st.selectbox("🎯 请选择要折叠的候选片段:", mol_options)
        
        if st.button("🏗️ 启动 ESMFold 极速建模", type="primary"):
            with st.spinner("🚀 正在将序列发送至云端计算集群，预测原子级构象 (通常需要 5-15 秒)..."):
                # 匹配选中的序列
                target_seq = ""
                for r in all_results:
                    if f"{r['序列名称 (FASTA ID)']} ({r['结构域']})" == selected_mol_name:
                        target_seq = r['完整序列']
                        break
                
                # 获取结构并渲染
                pdb_data = fetch_esm_fold_pdb(target_seq)
                if pdb_data:
                    st.success(f"✅ 建模成功！({selected_mol_name})")
                    col1, col2 = st.columns([2, 1])
                    with col1:
                        # 在网页左侧渲染 3D 图
                        render_3d_structure(pdb_data)
                    with col2:
                        st.markdown("#### 🖥️ 桌面端联动")
                        st.write("将该 PDB 文件导入 PyMOL、ChimeraX 或 Discovery Studio 以进行高精度表面电荷分析或分子对接。")
                        st.download_button(
                            label="📥 下载 .pdb 结构文件",
                            data=pdb_data,
                            file_name=f"{selected_mol_name.replace(' ', '_')}.pdb",
                            mime="protein/x-pdb",
                            type="primary"
                        )
                else:
                    st.error("❌ 建模失败：服务器拥挤或序列格式不被支持。")

# 👆 3D 实验室前端代码结束 👆

# ==========================================
# 5. 后勤补给：Excel 转 FASTA 清洗器 (丢失补回模块)
# ==========================================
st.markdown("---")
st.title("🗂️ 内部数据预处理：Excel 智能转 FASTA 清洗器")
st.info("💡 如果你们自家的测序或研发数据是 Excel 表格，请先拖入此处提取为 FASTA，然后再粘贴到页面最上方进行深度解析。")

uploaded_excel = st.file_uploader("📤 请上传包含序列的 Excel 或 CSV 文件", type=['xlsx', 'xls', 'csv'], key="excel_uploader")
if uploaded_excel is not None:
    if st.button("🔄 一键清洗并转换为 FASTA", type="primary"):
        try:
            if uploaded_excel.name.endswith('.csv'): df_raw = pd.read_csv(uploaded_excel, header=None)
            else: df_raw = pd.read_excel(uploaded_excel, header=None)
            
            fasta_lines = []
            invalid_names = ['nan', 'None', '', 'Protein name', 'Humanized variants', 'Heavy chain', 'Light chain', 'Final leading sequences', 'Humanized sequence:']
            for index, row in df_raw.iterrows():
                name = str(row.iloc[0]).strip() if len(row) > 0 else ""
                hc = str(row.iloc[1]).strip() if len(row) > 1 else ""
                lc = str(row.iloc[2]).strip() if len(row) > 2 else ""
                
                if name in invalid_names or "sequence" in name.lower() or name.startswith('Unnamed:'): continue
                if hc and hc not in ['nan', 'None'] and len(hc) > 10: fasta_lines.append(f">{name}_Heavy_Chain\n{hc.upper()}")
                if lc and lc not in ['nan', 'None'] and len(lc) > 10: fasta_lines.append(f">{name}_Light_Chain\n{lc.upper()}")
            
            fasta_content = '\n'.join(fasta_lines)
            if fasta_content:
                st.success(f"✅ 清洗完成！提取出 {fasta_content.count('>')} 条链序列。你可以直接点击下方按钮下载，或者全选展开的文本复制到上方。")
                st.download_button(label="📥 下载提纯后的 .fasta 文件", data=fasta_content, file_name="Cleaned_Sequences.fasta", mime="text/plain", type="primary")
                with st.expander("📄 点击预览清洗后的 FASTA 代码"):
                    st.code(fasta_content[:1500] + "\n\n... (仅展示部分预览)", language="text")
            else:
                st.warning("⚠️ 未能从文件中提取出有效序列，请检查表格格式。")
        except Exception as e:
            st.error(f"❌ 读取错误: {e}")

# ==========================================
# 6. 后勤补给：WIPO 专利序列表智能提取器
# ==========================================
st.markdown("---")
st.title("🗂️ 竞争情报获取：WIPO 专利序列表提取器")
st.info("💡 将 WIPO 下载的老版 ST.25 (.txt) 或新版 ST.26 (.xml) 专利序列表拖入，系统自动滤除核酸，极速提取带有专利专属前缀的纯蛋白质序列 Excel 库。")

uploaded_patent = st.file_uploader("📤 请上传专利序列表文件 (.txt 或 .xml)", type=['txt', 'xml'], key="patent_uploader")

if uploaded_patent is not None:
    default_patent_id = uploaded_patent.name.rsplit('.', 1)[0]
    patent_id = st.text_input("📝 确认专利号或项目名 (自动作为所有序列名称的前缀)", value=default_patent_id)

    if st.button("🔨 一键提取蛋白质序列", type="primary"):
        parsed_data = []
        file_extension = uploaded_patent.name.split('.')[-1].lower()
        
        try:
            # 引擎 A：ST.26 最新 XML 破译引擎
            if file_extension == 'xml':
                st.info("🔍 检测到 ST.26 XML 格式，启动深层解析...")
                tree = ET.parse(uploaded_patent)
                root = tree.getroot()
                
                for elem in root.iter():
                    if '}' in elem.tag: elem.tag = elem.tag.split('}', 1)[1]
                        
                for seq_data in root.findall('.//SequenceData'):
                    seq_id = seq_data.get('sequenceIDNumber', 'Unknown')
                    insd_seq = seq_data.find('.//INSDSeq')
                    if insd_seq is not None:
                        moltype = insd_seq.findtext('INSDSeq_moltype')
                        sequence = insd_seq.findtext('INSDSeq_sequence')
                        if moltype and moltype.upper() == 'AA' and sequence:
                            clean_seq = re.sub(r'[^A-Za-z]', '', sequence).upper()
                            if len(clean_seq) > 20: 
                                parsed_data.append({
                                    "Protein_Name": f"{patent_id}_SEQ_{seq_id}",
                                    "Heavy_Chain": clean_seq, "Light_Chain": "" 
                                })

            # 引擎 B：ST.25 经典 TXT 破译引擎
            elif file_extension == 'txt':
                st.info("🔍 检测到 ST.25 TXT 格式，启动多肽捕获...")
                content = uploaded_patent.getvalue().decode("utf-8", errors="ignore")
                blocks = re.split(r'<210>', content)[1:]
                
                for block in blocks:
                    seq_id_match = re.search(r'^\s*(\d+)', block)
                    if not seq_id_match: continue
                    seq_id = seq_id_match.group(1)
                    if '<212>  PRT' not in block and '<212> PRT' not in block: continue 
                        
                    seq_match = re.search(r'<400>.*?\n(.*?)<210>', block + '<210>', re.DOTALL)
                    if seq_match:
                        clean_seq = re.sub(r'[^A-Za-z]', '', seq_match.group(1)).upper()
                        if len(clean_seq) > 20:
                            parsed_data.append({
                                "Protein_Name": f"{patent_id}_SEQ_{seq_id}",
                                "Heavy_Chain": clean_seq, "Light_Chain": "" 
                            })

            # 渲染与导出
            if parsed_data:
                df_patent = pd.DataFrame(parsed_data)
                st.success(f"✅ 破译成功！共抓取到 **{len(parsed_data)}** 条合规蛋白质序列。")
                
                with st.expander("👀 预览前 10 条抓取序列", expanded=True):
                    st.dataframe(df_patent.head(10), use_container_width=True)
                
                buffer_patent = io.BytesIO()
                with pd.ExcelWriter(buffer_patent, engine='openpyxl') as writer:
                    df_patent.to_excel(writer, index=False, header=False)
                
                out_filename = f"DB_{patent_id}_Extracted.xlsx"
                st.download_button(
                    label=f"📥 下载 {patent_id} 专利弹药库 (.xlsx)",
                    data=buffer_patent.getvalue(),
                    file_name=out_filename,
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    type="primary"
                )
                st.caption("💡 提示：下载此表后，你可以直接将它拖入页面最下方的【FTO 专利侵权雷达】模块用作比对参考。")
            else:
                st.warning("⚠️ 扫描完成，但未在文件中发现合规蛋白质序列 (可能是纯核酸专利)。")

        except Exception as e:
            st.error(f"❌ 解析文件时发生错误: {e}")

# ==========================================
# 7. 终极防御：私有 FTO (自由实施权) 专利侵权雷达
# ==========================================
st.markdown("---")
st.title("🛡️ 独立模块：私有 FTO 专利侵权排查雷达")
st.info("💡 IP 防御阵地：上传你们整理的或刚才上方生成的『竞争对手核心专利序列库』。系统将自动与你页面最顶端输入的候选序列进行交叉同源性比对，提前锁定侵权地雷！")

col1, col2 = st.columns([1, 1])
with col1:
    db_file = st.file_uploader("📂 第一步：上传参考专利库 (.xlsx / .csv)", type=['xlsx', 'xls', 'csv'], key="fto_db")
    st.caption("提示：无表头，第1列为专利/序列名，第2列重链，第3列轻链。")
with col2:
    st.markdown("**第二步：配置侵权判定阈值**")
    risk_threshold = st.slider("🚨 红色高危预警阈值 (%)", min_value=70, max_value=100, value=85, step=1)
    safe_threshold = st.slider("⚠️ 黄色相似预警阈值 (%)", min_value=50, max_value=100, value=75, step=1)

if st.button("⚔️ 启动全库 FTO 侵权扫描", type="primary"):
    if not db_file:
        st.error("❌ 请先上传私有专利序列库！")
    elif not raw_input or ">" not in raw_input:
        st.error("❌ 请先在页面最上方的输入框中粘贴需要评估的自家候选序列 (FASTA)！")
    else:
        try:
            if db_file.name.endswith('.csv'): df_db = pd.read_csv(db_file, header=None)
            else: df_db = pd.read_excel(db_file, header=None)
                
            db_sequences = {}
            for index, row in df_db.iterrows():
                name = str(row.iloc[0]).strip()
                hc = str(row.iloc[1]).strip() if len(row) > 1 else ""
                lc = str(row.iloc[2]).strip() if len(row) > 2 else ""
                if name in ['nan', 'None', ''] or "sequence" in name.lower() or name.startswith('Unnamed:'): continue
                if hc and hc not in ['nan', 'None'] and len(hc) > 20: db_sequences[f"{name}_HC"] = hc.upper()
                if lc and lc not in ['nan', 'None'] and len(lc) > 20: db_sequences[f"{name}_LC"] = lc.upper()

            st.success(f"✅ 成功加载专利数据库：包含 {len(db_sequences)} 条参考链。")
            query_dict = parse_fasta(raw_input) 
            
            fto_results = []
            progress_bar = st.progress(0)
            total_comparisons = len(query_dict) * len(db_sequences)
            current_comp = 0
            
            for q_name, q_seq in query_dict.items():
                for db_name, db_seq in db_sequences.items():
                    matcher = difflib.SequenceMatcher(None, q_seq, db_seq)
                    homology_ratio = round(matcher.ratio() * 100, 2)
                    if homology_ratio >= safe_threshold:
                        risk_level = "🚨 极高危 (极可能侵权)" if homology_ratio >= risk_threshold else "⚠️ 中度风险 (结构相似)"
                        fto_results.append({
                            "我的候选分子": q_name, "匹配到的专利分子": db_name,
                            "序列同源性 (%)": homology_ratio, "风险定性": risk_level,
                            "候选序列长度": len(q_seq), "专利序列长度": len(db_seq)
                        })
                    current_comp += 1
                    if current_comp % 100 == 0 or current_comp == total_comparisons:
                        progress_bar.progress(min(current_comp / total_comparisons, 1.0))

            if fto_results:
                df_fto = pd.DataFrame(fto_results).sort_values(by="序列同源性 (%)", ascending=False)
                st.markdown("### 📑 FTO 侵权风险排查报告")
                def highlight_risk(val):
                    color = '#ff4b4b' if '高危' in str(val) else '#ffa726' if '中度' in str(val) else ''
                    return f'background-color: {color}'
                st.dataframe(df_fto.style.map(highlight_risk, subset=['风险定性']), use_container_width=True)
                
                buffer_fto = io.BytesIO()
                with pd.ExcelWriter(buffer_fto, engine='openpyxl') as writer:
                    df_fto.to_excel(writer, index=False, sheet_name='FTO_Risk_Report')
                st.download_button("📥 导出 FTO 排查报告 (.xlsx)", data=buffer_fto.getvalue(), file_name="FTO_Infringement_Report.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="primary")
            else:
                st.success("🎉 恭喜！你的所有候选序列与所上传的专利库同源性均低于安全阈值，FTO 风险极低，具有极好的自由实施空间！")
        except Exception as e:
            st.error(f"❌ 运行 FTO 雷达时发生错误: {e}")
