from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import streamlit as st
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor

Entrez.email = "your_email@example.com"

# DOIによる要約取得関数をキャッシュ
@lru_cache(maxsize=100)  # キャッシュサイズの指定
def fetch_summary_from_pubmed(doi):
    try:
        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={doi}"
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

def search_and_fetch_pubmed_articles(query, max_results=5):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, sort="most recent", retmax=max_results)
        record = Entrez.read(handle)
        id_list = record["IdList"]

        # 並列処理を導入
        with ThreadPoolExecutor(max_workers=max_results) as executor:
            articles = []
            futures = []
            for pubmed_id in id_list:
                summary_handle = Entrez.esummary(db="pubmed", id=pubmed_id)
                summary_record = Entrez.read(summary_handle)[0]
                
                title = summary_record.get("Title", "タイトルがありません")
                doi = summary_record.get("DOI", "DOIがありません")
                authors = ', '.join(summary_record.get("AuthorList", []))
                
                # 並列で要約取得を実行
                future = executor.submit(fetch_summary_from_pubmed, doi)
                futures.append((future, title, authors, doi))
            
            # 取得結果を順に処理
            for future, title, authors, doi in futures:
                abstract = future.result()  # 要約取得結果
                article = {
                    "title": title,
                    "authors": authors,
                    "abstract": abstract,
                    "doi": doi
                }
                articles.append(article)
        
        return articles
    except Exception as e:
        return []

# Streamlit UIの部分は変わりません
# Streamlitアプリのタイトルを設定
st.title("PubMed論文検索")

# 検索クエリを入力するテキストボックスを表示
query = st.text_input("検索クエリを入力してください")

# 検索ボタンがクリックされたら、PubMed APIを使用して論文を検索し、要約を取得して表示
if st.button("検索"):
    if query:
        articles = search_and_fetch_pubmed_articles(query, max_results=5)
        if articles:
            st.write(f"検索結果: {len(articles)} 件")
            for article in articles:
                st.write("-------")
                st.write(f"タイトル: {article['title']}")
                st.write(f"著者: {article['authors']}")
                st.write(f"要約: {article['abstract']}")
                st.write(f"DOI: {article['doi']}")
        else:
            st.warning("検索結果が見つかりませんでした。別のクエリをお試しください。")
    else:
        st.warning("検索クエリを入力してください。")