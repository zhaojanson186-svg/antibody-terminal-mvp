import streamlit as st
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from abnumber import Chain
import pandas as pd

st.set_page_config(page_title="抗体开发终端 | 智能分析", page_icon="🧬", layout="wide")

st.title("🧬 抗体序列智能分析中台 (防呆增强版)")
st.markdown("支持随意粘贴纯序列或 FASTA 格式，自动清洗非氨基酸字符，并批量提取 CDR/FR 区。")

with st.sidebar:
    st.header("⚙️ 设置")
    numbering_scheme = st.selectbox("选择编号规则", ("imgt", "chothia", "kabat", "martin"))

sequence_input = st.text_area(
    "📥 请在此处粘贴序列 (支持带 > 的 FASTA 格式或多条序列):", 
    height=250, 
    placeholder=">H16-7.8_Light_Chain (SEQ ID NO: 8)\nMLPSQLIGFLLL... (随便粘贴，程序会自动清洗)"
)

# --- 核心：智能文本清洗与解析函数 ---
def parse_dirty_input(text):
    records = []
    # 如果检测到 '>'，按 FASTA 逻辑处理多条
    if '>' in text:
        blocks = text.split('>')
        for block in blocks:
            if not block.strip(): continue
            lines = block.strip().split('\n')
            header = lines[0].strip()
            # 拼接剩余行，并用正则去除非字母字符（剔除数字、空格、标点）
            raw_seq = ''.join(lines[1:])
            clean_seq = re.sub(r'[^A-Za-z]', '', raw_seq).upper()
            if clean_seq:
                records.append({"id": header[:40], "seq": clean_seq})
    else:
        # 如果没有 '>'，当做单条序列处理，强力清洗
        clean_seq = re.sub(r'[^A-Za-z]', '', text).upper()
        if clean_seq:
            records.append({"id": "Input_Sequence", "seq": clean_seq})
    return records

if st.button("🚀 智能解析", type="primary"):
    if sequence_input.strip():
        records = parse_dirty_input(sequence_input)
        
        if not records:
            st.error("❌ 没有在输入中找到有效的氨基酸序列 (A-Z)。")
        else:
            st.success(f"✅ 成功清洗并识别出 {len(records)} 条序列！正在扫描可变区...")
            
            results = []
            progress_bar = st.progress(0)
            
            for i, rec in enumerate(records):
                seq_id = rec["id"]
                seq_str = rec["seq"]
                
                # 计算理化性质
                try:
                    pa = ProteinAnalysis(seq_str)
                    mw = round(pa.molecular_weight(), 2)
                    pi = round(pa.isoelectric_point(), 2)
                except:
                    mw, pi = 0.0, 0.0
                    
                # 识别抗体区域
                try:
                    chain = Chain(seq_str, scheme=numbering_scheme)
                    results.append({
                        "序列 ID": seq_id,
                        "链": chain.chain_type,
                        "MW (Da)": mw,
                        "pI": pi,
                        "FR1": chain.fr1_seq, "CDR1": chain.cdr1_seq,
                        "FR2": chain.fr2_seq, "CDR2": chain.cdr2_seq,
                        "FR3": chain.fr3_seq, "CDR3": chain.cdr3_seq,
                        "FR4": chain.fr4_seq,
                        "状态": "✅ 成功"
                    })
                except Exception as e:
                    results.append({
                        "序列 ID": seq_id, "链": "-", "MW (Da)": mw, "pI": pi,
                        "FR1": "-", "CDR1": "-", "FR2": "-", "CDR2": "-", 
                        "FR3": "-", "CDR3": "-", "FR4": "-",
                        "状态": "❌ 未找到完整V区"
                    })
                
                progress_bar.progress((i + 1) / len(records))

            # 渲染结果表格
            df = pd.DataFrame(results)
            
            # 美化表格：给 CDR 区域上色
            def style_df(x):
                df_style = pd.DataFrame('', index=x.index, columns=x.columns)
                for col in ['CDR1', 'CDR2', 'CDR3']:
                    df_style[col] = 'background-color: #e6f2ff; color: #004085'
                fail_mask = x['状态'].str.contains('❌')
                df_style.loc[fail_mask, '状态'] = 'color: red; font-weight: bold'
                return df_style

            st.dataframe(df.style.apply(style_df, axis=None), use_container_width=True)
            
            # 导出 CSV
            st.download_button(
                label="📥 导出分析报告 (CSV)",
                data=df.to_csv(index=False).encode('utf-8-sig'),
                file_name='Antibody_Report.csv',
                mime='text/csv',
            )
    else:
        st.warning("⚠️ 内容为空。")
