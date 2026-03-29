import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import io
from abnumber import Chain
import pandas as pd

st.set_page_config(page_title="抗体开发终端 | 多序列分析", page_icon="🧬", layout="wide")

st.title("🧬 抗体序列智能分析中台 (V2 多序列版)")
st.markdown("""
**支持功能：**
1. **多序列批量处理**：支持标准 FASTA 格式输入。
2. **全长序列智能截取**：自动忽略信号肽和恒定区 (Fc)，精准提取可变区 (V-domain)。
3. **理化性质计算**：计算全长序列的分子量和等电点。
""")

with st.sidebar:
    st.header("⚙️ 设置")
    numbering_scheme = st.selectbox(
        "选择编号规则",
        ("imgt", "chothia", "kabat", "martin")
    )

# 提示用户使用 FASTA 格式
sequence_input = st.text_area(
    "请输入氨基酸序列 (支持单条纯文本，或多条 FASTA 格式):", 
    height=250, 
    placeholder=""">Seq_1
MKHLWFFLLLVAAPRWVLSQVQLQESGPGLVKPSQTLSLTCTVSGGSISSGGYYWSWIRQ...
>Seq_2
DIQMTQSPSSLSASVGDRVTITCRASQGIRNDLGWYQQKPGKAPKRLIYAASSLQSGVPS..."""
)

if st.button("🚀 批量提取与分析", type="primary"):
    if sequence_input.strip():
        # 1. 智能解析输入内容
        records = []
        if ">" in sequence_input:
            # 识别为 FASTA 格式
            with io.StringIO(sequence_input) as f:
                for record in SeqIO.parse(f, "fasta"):
                    records.append({"id": record.id, "seq": str(record.seq).upper()})
        else:
            # 识别为单条纯文本
            clean_seq = "".join(sequence_input.split()).upper()
            records.append({"id": "Input_Seq", "seq": clean_seq})
            
        st.success(f"✅ 成功读取 {len(records)} 条序列，正在调用底层 HMMER 引擎进行扫描...")
        
        # 2. 批量处理
        results = []
        progress_bar = st.progress(0)
        
        for i, rec in enumerate(records):
            seq_id = rec["id"]
            seq_str = rec["seq"]
            
            # 理化性质（基于全长输入）
            try:
                pa = ProteinAnalysis(seq_str)
                mw = pa.molecular_weight()
                pi = pa.isoelectric_point()
            except:
                mw, pi = 0.0, 0.0
                
            # 抗体区域识别（自动在全长中寻找 V 区）
            try:
                chain = Chain(seq_str, scheme=numbering_scheme)
                results.append({
                    "序列 ID": seq_id,
                    "链类型": chain.chain_type,
                    "全长 MW (Da)": round(mw, 2),
                    "全长 pI": round(pi, 2),
                    "FR1": chain.fr1_seq,
                    "CDR1": chain.cdr1_seq,
                    "FR2": chain.fr2_seq,
                    "CDR2": chain.cdr2_seq,
                    "FR3": chain.fr3_seq,
                    "CDR3": chain.cdr3_seq,
                    "FR4": chain.fr4_seq,
                    "解析状态": "✅ 成功"
                })
            except Exception as e:
                # 记录解析失败的序列，但不中断整个程序
                results.append({
                    "序列 ID": seq_id,
                    "链类型": "-",
                    "全长 MW (Da)": round(mw, 2),
                    "全长 pI": round(pi, 2),
                    "FR1": "-", "CDR1": "-", "FR2": "-", "CDR2": "-", "FR3": "-", "CDR3": "-", "FR4": "-",
                    "解析状态": f"❌ 失败 (未找到完整V区)"
                })
            
            progress_bar.progress((i + 1) / len(records))

        # 3. 结果展示
        if results:
            df = pd.DataFrame(results)
            
            # 统计汇总面板
            success_count = len(df[df['解析状态'] == '✅ 成功'])
            st.markdown("### 📊 分析汇总")
            st.write(f"共处理 **{len(records)}** 条序列，成功识别出 **{success_count}** 个抗体可变区。")
            
            # 高亮显示 CDR 区域的函数
            def style_dataframe(x):
                # 仅对成功的行，将 CDR 列标蓝
                df_style = pd.DataFrame('', index=x.index, columns=x.columns)
                for col in ['CDR1', 'CDR2', 'CDR3']:
                    df_style[col] = 'background-color: #e6f2ff; color: #004085'
                # 失败状态标红
                fail_mask = x['解析状态'].str.contains('❌')
                df_style.loc[fail_mask, '解析状态'] = 'color: red; font-weight: bold'
                return df_style

            st.dataframe(df.style.apply(style_dataframe, axis=None), use_container_width=True)
            
            # 提供 CSV 下载
            csv = df.to_csv(index=False).encode('utf-8-sig')
            st.download_button(
                label="📥 导出结果为 CSV (支持 Excel 打开)",
                data=csv,
                file_name='antibody_analysis_results.csv',
                mime='text/csv',
            )

    else:
        st.warning("⚠️ 请先输入序列！")
