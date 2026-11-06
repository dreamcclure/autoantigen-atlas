"""
Phase 8 — Autoantigen Atlas web app (hosted chat + discovery + 3D structures).
=============================================================================
A browser interface over the existing engine, in three tabs:
  • 💬 Ask the Atlas  — chat with the cited agent (wraps src/agent.py)
  • 🔬 Discovery      — ranked candidate autoantigen hypotheses (reads src/discover.py output)
  • 🧊 3D Structures  — the AlphaFold viewer embedded in-app (figures/*.html)
No engine is re-implemented here: the agent, the discovery runs, and the viewers all come
from code and data produced elsewhere in the project. This file is pure presentation.

Run locally:   streamlit run app.py     (needs ANTHROPIC_API_KEY in .env)
Deploy:        Streamlit Community Cloud — put the key in Secrets, and set a hard
               monthly spend cap in the Anthropic console (owner-funded demo).
"""

import json
import os
import sys
import time

import streamlit as st
import streamlit.components.v1 as components

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
PROCESSED = os.path.join(HERE, "data", "processed")
FIGURES = os.path.join(HERE, "figures")

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


# ============================================================ agent (tab 1)
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


# ============================================================ discovery (tab 2)
@st.cache_data(show_spinner=False)
def load_discovery(disease):
    """Load the precomputed partner-arm and (if present) literature-arm discovery runs."""
    def _read(name):
        try:
            with open(os.path.join(PROCESSED, name)) as f:
                return json.load(f)
        except FileNotFoundError:
            return None
    return (_read("discovery_{}.json".format(disease)),
            _read("discovery_{}_literature.json".format(disease)))


# evidence tiers the agent assigns, with a display style (emoji + Streamlit color)
TIER = {
    "ESTABLISHED": ("🟢", "green", "Established"),
    "EMERGING":    ("🟡", "orange", "Emerging"),
    "NONE-FOUND":  ("⚪", "gray", "No published evidence"),
}
_TIER_ORDER = {"ESTABLISHED": 0, "EMERGING": 1, "NONE-FOUND": 2}


def _na_label(value):
    return {True: "yes", False: "no", None: "—"}.get(value, "—")


def _render_assessments(assessments):
    """Show the agent's cited verdict per candidate, strongest evidence first."""
    for a in sorted(assessments, key=lambda x: _TIER_ORDER.get(x["evidence_tier"], 9)):
        emoji, color, label = TIER.get(a["evidence_tier"], ("", "gray", a["evidence_tier"]))
        st.markdown("**{}** (`{}`) — {} :{}[{}]".format(
            a["gene"], a["accession"], emoji, color, label))
        with st.expander("Agent's cited assessment"):
            st.markdown(a["transcript"])


def render_discovery(disease):
    partner, literature = load_discovery(disease)
    st.caption("Investigating **{}**.".format(DISEASES[disease]["name"]))
    st.info("**Discovery mode flips the question** — instead of *“is protein X a known "
            "autoantigen?”* it asks *“which proteins, **not** already established, might be?”* "
            "Everything below is a **ranked, evidence-backed hypothesis for expert review — "
            "never a confirmed discovery.**")

    if not partner and not literature:
        st.warning("No discovery run is available for this disease yet.")
        return

    # ---- Partner arm: epitope-spreading neighbors of known antigens ----
    if partner and partner.get("candidates"):
        st.markdown("#### 🧩 Partner arm — epitope-spreading neighbors")
        with st.expander("How this arm works"):
            st.caption(partner["_metadata"]["method"])
        assessed = {a["accession"]: a for a in partner.get("agent_assessments", [])}
        rows = []
        for c in partner["candidates"][:15]:
            rows.append({
                "Gene": c["gene"],
                "Linked antigens": c.get("partner_support"),
                "Specificity": round(c.get("specificity", 0), 3),
                "Nuclear/NA": _na_label(c.get("na_associated")),
                "Connects to": ", ".join(c.get("connects_to", [])[:4]),
            })
        st.dataframe(rows, width="stretch", hide_index=True)
        st.caption("*Specificity* = fraction of this protein's partners that are known antigens "
                   "(controls for promiscuous interaction hubs). Showing top "
                   "{} of {} candidates.".format(min(15, len(partner["candidates"])),
                                                 len(partner["candidates"])))
        if partner.get("agent_assessments"):
            st.markdown("**Agent assessment of the top candidates:**")
            _render_assessments(partner["agent_assessments"])

    # ---- Literature arm: autoantibody-context co-mentions ----
    if literature and literature.get("candidates"):
        st.divider()
        st.markdown("#### 📚 Literature arm — autoantibody-context co-mentions")
        with st.expander("How this arm works"):
            st.caption(literature["_metadata"]["method"])
        assessed = {a["accession"]: a for a in literature.get("agent_assessments", [])}
        rows = []
        for c in literature["candidates"]:
            a = assessed.get(c["accession"])
            tier = TIER.get(a["evidence_tier"])[2] if a else "—"
            rows.append({
                "Gene": c["gene"],
                "Accession": c["accession"],
                "Papers": c.get("lit_support"),
                "Agent verdict": tier,
            })
        st.dataframe(rows, width="stretch", hide_index=True)
        if literature.get("agent_assessments"):
            st.markdown("**Agent assessment of the candidates:**")
            _render_assessments(literature["agent_assessments"])


# ============================================================ 3D structures (tab 3)
VIEWER = {"lupus": "structure_viewer.html",
          "rheumatoid_arthritis": "rheumatoid_arthritis_viewer.html"}


def render_structures(disease):
    st.caption("Every structure is an **AlphaFold** prediction, colored by per-residue confidence "
               "(**pLDDT**): blue = very high, orange/red = low (often flexible / disordered "
               "regions). Drag to rotate.")
    name = VIEWER.get(disease)
    if not name:
        st.info("No 3D viewer is available for this disease yet.")
        return
    # Embed the viewer by its hosted URL rather than inlining the file: the viewer is a full
    # 100vh canvas app, and a real iframe gives it a proper viewport height (inlining it left
    # the canvas measuring 0 height). GitHub Pages sends no X-Frame-Options, so framing is OK.
    url = "{}/{}".format(PAGES, name)
    components.iframe(url, height=820, scrolling=True)
    st.caption("Trouble loading the embed? [Open the full-screen viewer in a new tab]({}).".format(url))


# ============================================================ sidebar
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
    st.markdown("🔬 **Open 3D viewers in a new tab**  \n"
                "[Lupus]({}/structure_viewer.html) · "
                "[RA]({}/rheumatoid_arthritis_viewer.html)".format(PAGES, PAGES))
    st.markdown("💻 [Source on GitHub]({})".format(REPO))
    st.caption("Research aide for scientists — not medical advice.")

# ============================================================ main
st.title("🧬 Autoantigen Atlas")

# reset the conversation when the disease changes
if st.session_state.get("disease") != disease:
    st.session_state.disease = disease
    st.session_state.messages = []

tab_chat, tab_discovery, tab_structures = st.tabs(
    ["💬 Ask the Atlas", "🔬 Discovery", "🧊 3D Structures"])

# ---------------------------------------------------------------- tab 1: chat
with tab_chat:
    st.caption("Investigating **{}**. Ask whether a protein is an autoantigen, what its structure "
               "looks like, or what the literature says — every answer comes with citations."
               .format(DISEASES[disease]["name"]))

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
            if col.button(ex, width="stretch"):
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
                st.warning("You've reached this demo's hourly question limit. Please try again later "
                           "— or run it locally with your own API key (see the GitHub repo).")
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

# ---------------------------------------------------------------- tab 2: discovery
with tab_discovery:
    render_discovery(disease)

# ---------------------------------------------------------------- tab 3: 3D structures
with tab_structures:
    render_structures(disease)
