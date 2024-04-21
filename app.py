import requests
from bs4 import BeautifulSoup
import streamlit as st
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor

# 直接DeepL APIキーを入力（セキュリティのため、このような直接入力は避ける方がよい）
DEEPL_API_KEY = 'dbe78709-8805-4d30-be62-eb4acc1b5392:fx'

if not DEEPL_API_KEY:
    raise ValueError("DeepL APIキーが設定されていません")

Entrez.email = "your_email@example.com"

# DeepL APIを使ってテキストを翻訳する関数
def translate_to_japanese(text):
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

# PubMedからDOIを取得して要約を抽出する関数
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
        return f"エラーが発生しました: {e}"

# PubMed論文を検索して情報を取得する関数
def search_and_fetch_pubmed_articles(query, max_results=5):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, sort="most recent", retmax=max_results)
        record = Entrez.read(handle)
        id_list = record["IdList"]

        articles = []
        with ThreadPoolExecutor(max_workers=max_results) as executor:
            futures = []
            for pubmed_id in id_list:
                future = executor.submit(fetch_abstract_from_pubmed, pubmed_id)
                futures.append((future, pubmed_id))
            
            for future, pubmed_id in futures:
                abstract = future.result()
                translated_abstract = translate_to_japanese(abstract)
                
                article = {
                    "title": f"PubMed Article {pubmed_id}",
                    "abstract": translated_abstract,
                    "pubmed_id": pubmed_id
                }
                articles.append(article)
        
        return articles
    except Exception as e:
        return f"エラーが発生しました: {e}"

# Streamlitアプリケーション
st.title("PubMed論文検索")

query = st.text_input("検索クエリを入力してください")

if st.button("検索"):
    if query:
        articles = search_and_fetch_pubmed_articles(query, max_results=5)
        if articles:
            st.write(f"検索結果: {len(articles)} 件")
            for article in articles:
                st.write("-------")
                st.write(f"タイトル: {article['title']}")
                st.write(f"要約: {article['abstract']}")
                st.write(f"PubMed ID: {article['pubmed_id']}")
        else:
            st.warning("検索結果が見つかりませんでした。別のクエリをお試しください。")
    else:
        st.warning("検索クエリを入力してください。")
