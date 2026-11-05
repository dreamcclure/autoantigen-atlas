"""
Phase 8 — Autoantigen Atlas web app (hosted chat).
==================================================
A browser chat interface over the existing agent: pick a disease, ask about a protein,
read a CITED answer. Wraps src/agent.py — no re-implementation of the engine.

Run locally:   streamlit run app.py     (needs ANTHROPIC_API_KEY in .env)
Deploy:        Streamlit Community Cloud — put the key in Secrets, and set a hard
               monthly spend cap in the Anthropic console (owner-funded demo).
"""

import os
import sys
import time

import streamlit as st

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))

# Secrets: local uses .env; Streamlit Cloud injects st.secrets → bridge it to the env
# so ChatAnthropic (which reads ANTHROPIC_API_KEY from the environment) finds the key.
try:
    from dotenv import load_dotenv
    load_dotenv(os.path.join(HERE, ".env"))
except Exception:
    pass
try:
    if "ANTHROPIC_API_KEY" in st.secrets:
        os.environ["ANTHROPIC_API_KEY"] = st.secrets["ANTHROPIC_API_KEY"]
except Exception:
    pass

from diseases import DISEASES

# --- cost guardrail: per-session rate limit (the hard ceiling is the console spend cap) ---
RATE_LIMIT = 10        # questions...
RATE_WINDOW = 3600     # ...per hour, per browser session

st.set_page_config(page_title="Autoantigen Atlas", page_icon="🧬", layout="centered")


@st.cache_resource(show_spinner=False)
def get_agent(disease):
    """Build the agent once per disease (cached across reruns / messages)."""
    from agent import build_agent
    return build_agent(disease)


def _text(content):
    """Normalize an assistant message's content (string or list-of-blocks) to text."""
    if isinstance(content, list):
        return "".join(b.get("text", "") for b in content
                       if isinstance(b, dict) and b.get("type") == "text")
    return content or ""


def ask_agent(agent, question):
    """One agent run. Returns (answer_text, [tool names used], error_or_None)."""
    try:
        result = agent.invoke({"messages": [{"role": "user", "content": question}]})
    except Exception as error:
        return None, [], str(error)
    messages = result["messages"]
    tools_used = []
    for m in messages:
        for call in getattr(m, "tool_calls", None) or []:
            tools_used.append(call.get("name", "?"))
    return _text(messages[-1].content), tools_used, None


def within_rate_limit():
    now = time.time()
    hits = [t for t in st.session_state.get("hits", []) if now - t < RATE_WINDOW]
    st.session_state.hits = hits
    return len(hits) < RATE_LIMIT


# ---------------------------------------------------------------- sidebar
TOOL_LABEL = {"lookup_protein": "🧬 UniProt", "get_structure": "🧊 AlphaFold",
              "search_literature": "🔎 Europe PMC"}
REPO = "https://github.com/dreamcclure/autoantigen-atlas"
PAGES = "https://dreamcclure.github.io/autoantigen-atlas/figures"

with st.sidebar:
    st.header("Autoantigen Atlas")
    st.caption("An AI research aide that finds the proteins the immune system attacks in "
               "autoimmune disease — and **cites every claim**.")
    disease = st.selectbox("Disease", list(DISEASES.keys()),
                           format_func=lambda k: DISEASES[k]["name"])
    st.divider()
    st.markdown(
        "**How it works** — a Claude Opus 4.8 agent with three free tools "
        "(UniProt, AlphaFold, Europe PMC) and five guardrails: it answers **only** from "
        "tool results, cites every source, and declines rather than guesses.")
    st.markdown("**Blind-tested:** lupus 96% recall · RA 81% · 100% precision & specificity on both.")
    st.divider()
    st.markdown("🔬 **3D structure viewers**  \n"
                "[Lupus]({}/structure_viewer.html) · "
                "[RA]({}/rheumatoid_arthritis_viewer.html)".format(PAGES, PAGES))
    st.markdown("💻 [Source on GitHub]({})".format(REPO))
    st.caption("Research aide for scientists — not medical advice.")

# ---------------------------------------------------------------- main
st.title("🧬 Ask the Atlas")
st.caption("Investigating **{}**. Ask whether a protein is an autoantigen, what its structure "
           "looks like, or what the literature says — every answer comes with citations."
           .format(DISEASES[disease]["name"]))

# reset the conversation when the disease changes
if st.session_state.get("disease") != disease:
    st.session_state.disease = disease
    st.session_state.messages = []

if not st.session_state.get("messages"):
    st.markdown("**Try asking:**")
    examples = {
        "lupus": ["Is Ro60 an autoantigen in lupus?",
                  "What does the literature say about anti-C1q antibodies?",
                  "How confident is the AlphaFold structure for La/SSB?"],
        "rheumatoid_arthritis": ["Is vimentin an autoantigen in rheumatoid arthritis?",
                                 "What's the evidence for anti-citrullinated protein antibodies?",
                                 "Tell me about alpha-enolase in RA."],
    }.get(disease, [])
    cols = st.columns(len(examples)) if examples else []
    for col, ex in zip(cols, examples):
        if col.button(ex, use_container_width=True):
            st.session_state.pending = ex

# replay history
for msg in st.session_state.get("messages", []):
    with st.chat_message(msg["role"]):
        st.markdown(msg["content"])
        if msg.get("tools"):
            st.caption("Consulted: " + " · ".join(TOOL_LABEL.get(t, t) for t in msg["tools"]))

prompt = st.chat_input("Ask about a protein…") or st.session_state.pop("pending", None)

if prompt:
    st.session_state.messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    if not within_rate_limit():
        with st.chat_message("assistant"):
            st.warning("You've reached this demo's hourly question limit. Please try again later — "
                       "or run it locally with your own API key (see the GitHub repo).")
    else:
        st.session_state.hits = st.session_state.get("hits", []) + [time.time()]
        with st.chat_message("assistant"):
            with st.status("Investigating — consulting UniProt, AlphaFold, and Europe PMC…",
                           expanded=False):
                answer, tools, error = ask_agent(get_agent(disease), prompt)
            if error:
                st.error("The agent hit a problem: {}. Please try again.".format(error))
            else:
                st.markdown(answer)
                if tools:
                    st.caption("Consulted: " + " · ".join(
                        dict.fromkeys(TOOL_LABEL.get(t, t) for t in tools)))
                st.session_state.messages.append(
                    {"role": "assistant", "content": answer, "tools": tools})
