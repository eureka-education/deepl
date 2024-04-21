import os
import requests
from bs4 import BeautifulSoup
import streamlit as st
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor
import re

# DeepL APIキーとEntrezメールアドレスはシークレットから取得
DEEPL_API_KEY = st.secrets["DEEPL_API_KEY"]
ENTREZ_EMAIL = st.secrets["ENTREZ_EMAIL"]

# Entrezでのリクエスト用メールアドレスを設定
Entrez.email = ENTREZ_EMAIL

# DeepL APIを使ってテキストを翻訳する関数
def translate_to_japanese(text):
    try:
        url = "https://api-free.deepl.com/v2/translate"
        payload = {
            "auth_key": DEEPL_API_KEY,
            "text": text,
            "target_lang": "JA"
        }
        response = requests.post(url, data=payload)

        if response.status_code == 200:
            translated_data = response.json()
            translated_text = translated_data["translations"][0]["text"]
            return translated_text
        else:
            return f"翻訳に失敗しました。ステータスコード: {response.status_code}, メッセージ: {response.text}"
    except Exception as e:
        return f"DeepL翻訳でエラーが発生しました: {e}"

# PubMedから要約を抽出する関数
def fetch_abstract_from_pubmed(pubmed_id):
    try:
        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
        response = requests.get(pubmed_url)
        if response.status_code == 200:
            soup = BeautifulSoup(response.content, "html.parser")
            abstract_div = soup.find("div", class_="abstract-content")
            if abstract_div:
                abstract_text = ' '.join(abstract_div.stripped_strings)
                return abstract_text
        return "要約がありません"
    except Exception as e:
        return f"要約取得でエラーが発生しました: {e}"

# PubMedで論文を検索して情報を取得する関数
def search_and_fetch_pubmed_articles(query, max_results=5):
    try:
        # PubMedで論文を検索
        handle = Entrez.esearch(db="pubmed", term=query, sort="most recent", retmax=max_results)
        record = Entrez.read(handle)
        id_list = record["IdList"]

        articles = []
        with ThreadPoolExecutor(max_workers=max_results) as executor:
            futures = []
            for pubmed_id in id_list:
                # 論文の概要情報を取得
                summary_handle = Entrez.esummary(db="pubmed", id=pubmed_id)
                summary_record = Entrez.read(summary_handle)[0]

                title = summary_record.get("Title", "タイトルがありません")
                authors = ', '.join(summary_record.get("AuthorList", ["著者がありません"]))

                # 要約を取得して翻訳
                future = executor.submit(fetch_abstract_from_pubmed, pubmed_id)
                futures.append((future, title, authors, pubmed_id))
            
            for future, title, authors, pubmed_id in futures:
                abstract = future.result()  # 原文の要約
                translated_abstract = translate_to_japanese(abstract)  # 翻訳された要約

                # 論文情報を作成
                article = {
                    "title": title,
                    "authors": authors,
                    "abstract": translated_abstract,
                    "pubmed_id": pubmed_id
                }
                articles.append(article)
        
        return articles
    except Exception as e:
        return f"PubMed検索でエラーが発生しました: {e}"

# Streamlitアプリケーション
st.title("PubMed論文検索")

# クエリのバリデーション
def is_valid_query(query):
    # アルファベット、数字、スペース、ハイフン、アンダースコアのみ許可
    return re.match(r'^[\w -]+$', query) is not None

query = st.text_input("検索クエリを入力してください")

if st.button("検索"):
    if is_valid_query(query):
        articles = search_and_fetch_pubmed_articles(query, max_results=5)
        if isinstance(articles, list):
            st.write(f"検索結果: {len(articles)} 件")
            for article in articles:
                st.write("-------")
                st.write(f"タイトル: {article['title']}")
                st.write(f"著者: {article['authors']}")
                st.write(f"要約: {article['abstract']}")
                st.write(f"PubMed ID: {article['pubmed_id']}")
        else:
            st.error(articles)  # エラーメッセージを表示
    else:
        st.warning("無効なクエリが入力されました。許可されている文字のみを使用してください。")

