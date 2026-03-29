import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from abnumber import Chain
import pandas as pd

# 页面配置
st.set_page_config(page_title="抗体开发终端 | 序列分析", page_icon="🧬", layout="wide")

st.title("🧬 抗体序列分析与注释工具 (MVP版)")
st.markdown("输入抗体重链或轻链氨基酸序列，自动计算理化性质并提取 CDR/FR 区域。")

# 侧边栏配置
with st.sidebar:
    st.header("⚙️ 设置")
    numbering_scheme = st.selectbox(
        "选择编号规则 (Numbering Scheme)",
        ("imgt", "chothia", "kabat", "martin")
    )
    st.markdown("---")
    st.markdown("💡 *提示: 推荐使用默认的 IMGT 规则进行 CDR 划分。*")

# 主界面输入
sequence_input = st.text_area("请输入单条氨基酸序列 (直接粘贴纯文本):", height=150, placeholder="例如: EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK...")

if st.button("🚀 开始分析", type="primary"):
    if sequence_input.strip():
        # 清理序列（去除空格和换行）
        clean_seq = "".join(sequence_input.split()).upper()
        
        col1, col2 = st.columns(2)
        
        # 模块 1: 理化性质计算
        with col1:
            st.subheader("📊 理化性质")
            try:
                analyzed_seq = ProteinAnalysis(clean_seq)
                mw = analyzed_seq.molecular_weight()
                pi = analyzed_seq.isoelectric_point()
                
                # 使用 metric 组件美观展示
                m1, m2 = st.columns(2)
                m1.metric("分子量 (Molecular Weight)", f"{mw:.2f} Da")
                m2.metric("等电点 (pI)", f"{pi:.2f}")
            except Exception as e:
                st.error(f"理化性质计算出错: 请确保输入的是标准氨基酸序列。")

        # 模块 2: 序列注释与 CDR 提取 (调用 ANARCI 底层)
        with col2:
            st.subheader(f"🏷️ 区域划分 (基于 {numbering_scheme.upper()})")
            try:
                # 实例化 Chain 对象，自动识别链类型并编号
                chain = Chain(clean_seq, scheme=numbering_scheme)
                
                st.info(f"**识别结果:** 这是属于 **{chain.chain_type}** 链的序列。")
                
                # 提取各个区域
                regions_data = {
                    "区域 (Region)": ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"],
                    "序列 (Sequence)": [
                        chain.fr1_seq, chain.cdr1_seq, chain.fr2_seq, 
                        chain.cdr2_seq, chain.fr3_seq, chain.cdr3_seq, chain.fr4_seq
                    ]
                }
                
                df_regions = pd.DataFrame(regions_data)
                
                # 在表格中高亮 CDR 区域以增强可读性
                def highlight_cdr(s):
                    if 'CDR' in s['区域 (Region)']:
                        return ['background-color: #e6f2ff'] * len(s)
                    return [''] * len(s)
                
                st.dataframe(df_regions.style.apply(highlight_cdr, axis=1), use_container_width=True)
                
            except Exception as e:
                st.error("序列解析失败。请检查序列是否完整（通常需要包含完整的 V 区）。")
                st.code(str(e))
    else:
        st.warning("⚠️ 请先输入序列！")
