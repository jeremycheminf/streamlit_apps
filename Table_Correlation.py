import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import io

st.title("Correlation and Pair Plot App")

st.write("Paste your tabular data below:")

data = st.text_area("paste data here", height=200, label_visibility="hidden")

if data:
    try:
        df = pd.read_csv(io.StringIO(data), sep="\t")
        st.write("Data Preview:")
        st.write(df.head())

        correlations = df.corr()
        st.write("Correlation Matrix:")
        st.write(correlations)

        st.write("Correlation Types:")
        correlation_type = st.selectbox(
            label="Select a correlation type",
            options=["Pearson", "Spearman"],
            label_visibility="hidden",
        )

        if correlation_type == "Pearson":
            correlations = df.corr(method="pearson")
        elif correlation_type == "Spearman":
            correlations = df.corr(method="spearman")

        st.write(correlations)

        st.write("Pair Plot:")
        fig = sns.pairplot(df)
        st.pyplot(fig)

    except Exception as e:
        st.write("Error parsing data: ", str(e))
