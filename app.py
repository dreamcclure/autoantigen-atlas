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
import discover   # reused for live candidate ranking + the assessment prompt/parser

# --- cost guardrail: per-session rate limit (the hard ceiling is the console spend cap) ---
RATE_LIMIT = 10        # questions...
RATE_WINDOW = 3600     # ...per hour, per browser session

st.set_page_config(page_title="Autoantigen Atlas", page_icon="🧬", layout="wide")

# A few polish touches the theme config can't express (gradient ground, heading sizes, borders).
# The fonts themselves are now set natively via theme.font / theme.headingFont (self-hosted woff2
# in static/) in .streamlit/config.toml — the old Google-Fonts @import here was silently ignored.
st.markdown("""
<style>
  .stApp { line-height: 1.6; }

  /* editorial title: stacked + uppercase, with ATLAS dropped below and offset (echoes the site) */
  .atlas-title { display: flex; flex-direction: column; align-items: flex-start;
    width: fit-content; text-transform: uppercase; line-height: 0.95; letter-spacing: -0.005em;
    font-weight: 600; font-size: 4.3rem !important; margin: 0.2rem 0 2.6rem !important; }
  /* ATLAS steps down and to the right, à la the DREA MCCLURE nameplate */
  .atlas-title .l2 { margin-left: 2.4em; letter-spacing: 0.1em; }

  /* larger, airier section headings (Fraunces face comes from theme.headingFont in config.toml) */
  h2 { font-size: 2.3rem !important; font-weight: 600 !important; letter-spacing: -0.01em;
       margin: 2.2rem 0 0.8rem !important; }
  h3 { font-size: 1.6rem !important; font-weight: 600 !important; margin: 1.6rem 0 0.6rem !important; }

  /* keep sidebar headings compact so the long title doesn't break mid-word in the narrow rail */
  section[data-testid="stSidebar"] h1,
  section[data-testid="stSidebar"] h2,
  section[data-testid="stSidebar"] h3 { font-size: 1.45rem !important; line-height: 1.2 !important;
       margin: 0.1rem 0 0.5rem !important; overflow-wrap: anywhere; }

  /* more breathing room around + between the main page's elements */
  .block-container { padding-top: 3.5rem !important; padding-bottom: 4rem !important; }
  section[data-testid="stMain"] [data-testid="stVerticalBlock"] { gap: 1.8rem !important; }
  /* Streamlit zeroes inner heading margins, so put real section breaks on the flex-child wrapper */
  section[data-testid="stMain"] [data-testid="stElementContainer"]:has(.sec-head) {
    margin-top: 3.6rem !important; }
  section[data-testid="stMain"] [data-testid="stElementContainer"]:has(.tryasking) {
    margin-top: 1rem !important; }
  section[data-testid="stMain"] [data-testid="stElementContainer"]:has([data-testid="stChatInput"]) {
    margin-top: 1.2rem !important; }

  /* the sidebar's Disease dropdown: a white box with dark text (the rail itself stays dark) */
  section[data-testid="stSidebar"] [data-baseweb="select"] > div {
    background-color: #FFFFFF !important; border-color: #FFFFFF !important; }
  section[data-testid="stSidebar"] [data-baseweb="select"] > div * { color: #1B1212 !important; }

  /* keep accessions/code monospaced */
  code, pre, kbd { font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }

  div[data-testid="stImage"] img { border-radius: 0.6rem; }

  /* main tabs: uppercase, evenly spread under the title, blue on hover */
  section[data-testid="stMain"] [data-baseweb="tab-list"] {
    width: 100%; justify-content: space-between; gap: 0; margin-bottom: 2.4rem; }
  section[data-testid="stMain"] button[data-baseweb="tab"],
  section[data-testid="stMain"] button[data-baseweb="tab"] p {
    text-transform: uppercase; letter-spacing: 0.1em; font-size: 0.95rem; font-weight: 600; }
  section[data-testid="stMain"] button[data-baseweb="tab"]:hover,
  section[data-testid="stMain"] button[data-baseweb="tab"]:hover p { color: #2563EB !important; }
  /* active tab stays licorice (blue only on hover), and no underline bar / divider */
  section[data-testid="stMain"] button[data-baseweb="tab"][aria-selected="true"],
  section[data-testid="stMain"] button[data-baseweb="tab"][aria-selected="true"] p { color: #1B1212 !important; }
  section[data-testid="stMain"] [data-baseweb="tab-highlight"],
  section[data-testid="stMain"] [data-baseweb="tab-border"] { background-color: transparent !important; }

  /* larger gap between the header row and the tab row (+ SCROLL) */
  section[data-testid="stMain"] [data-testid="stHorizontalBlock"]:has([data-baseweb="tab-list"]) {
    margin-top: 3rem !important; }

  /* space after the tab bar, before the intro line */
  section[data-testid="stMain"] [data-baseweb="tab-panel"] { padding-top: 1.4rem; }

  /* "Try asking" — Fraunces, uppercase, with space above it */
  .tryasking { font-family: 'Inter', sans-serif; text-transform: uppercase; color: #1B1212;
    font-size: 1rem; font-weight: 600; letter-spacing: 0.08em; margin: 2.6rem 0 0.5rem; }

  /* chat input: larger, a licorice bar with Fraunces text */
  [data-testid="stChatInput"] > div,
  [data-testid="stChatInput"] [data-baseweb="textarea"],
  [data-testid="stChatInput"] [data-baseweb="base-input"] {
    background: #1B1212 !important; border-color: #1B1212 !important; }
  [data-testid="stChatInput"] > div { border-radius: 0.7rem !important; }
  [data-testid="stChatInput"] [data-baseweb="base-input"] { padding: 0.55rem 0.4rem !important; }
  [data-testid="stChatInput"] textarea {
    color: #FFFFFF !important; font-family: 'Fraunces', Georgia, serif !important;
    font-size: 1.4rem !important; }
  [data-testid="stChatInput"] textarea::placeholder {
    color: rgba(255,255,255,0.8) !important; font-family: 'Fraunces', Georgia, serif !important; }
  [data-testid="stChatInput"] button svg { fill: #FDFDFD !important; color: #FDFDFD !important; }

  /* sidebar section headings — Fraunces, uppercase, on the dark rail */
  .side-head { font-family: 'Fraunces', Georgia, serif; text-transform: uppercase;
    font-size: 1.3rem; font-weight: 600; letter-spacing: 0.04em; color: #E6E9F2;
    margin: 0.3rem 0 0.6rem; }

  /* main Fraunces section titles (Autoimmune Disease, Protein Profile) */
  .sec-head { font-family: 'Fraunces', Georgia, serif; text-transform: uppercase;
    font-size: 1.55rem; font-weight: 600; letter-spacing: 0.04em; color: #1B1212;
    margin: 2.6rem 0 0.7rem; }

  /* the research-aide blurb to the right of the title */
  .hero-desc { font-family: 'Inter', sans-serif; font-weight: 300; font-size: 1.15rem;
    line-height: 1.5; color: #4a4443; max-width: 34ch; margin: 0; }
  .hero-desc b { font-weight: 600; color: #1B1212; }

  /* SCROLL cue — Fraunces, to the right of the tabs; jumps to the disease section */
  .scroll-link { display: inline-block; font-family: 'Fraunces', Georgia, serif;
    text-transform: uppercase; font-size: 1rem; font-weight: 600; letter-spacing: 0.12em;
    color: #1B1212 !important; text-decoration: none !important; white-space: nowrap; }
  .scroll-link:hover { color: #2563EB !important; }

  /* the disease selection bar in the main area — large */
  section[data-testid="stMain"] [data-baseweb="select"] > div {
    min-height: 3.3rem; font-size: 1.15rem; border-radius: 0.6rem; }

  /* Autoimmune Disease — a black "hero" band with white type (like the site hero) */
  .st-key-disease_band { background: #1B1212 !important; border-radius: 0.6rem !important;
    padding: 1.6rem 1.8rem 1.9rem !important; }
  .band-head { color: #FDFDFD; font-family: 'Fraunces', Georgia, serif; text-transform: uppercase;
    font-size: 1.7rem; font-weight: 600; letter-spacing: 0.04em; margin: 0 0 0.9rem; }
  .st-key-disease_band [data-baseweb="select"] > div {
    background: rgba(255,255,255,0.08) !important; border: 1px solid rgba(255,255,255,0.4) !important; }
  .st-key-disease_band [data-baseweb="select"] > div *,
  .st-key-disease_band [data-baseweb="select"] svg {
    color: #FDFDFD !important; fill: #FDFDFD !important; }


  /* air above the tab bar; wide layout left-aligned with a gap from the sidebar */
  section[data-testid="stMain"] [data-baseweb="tab-list"] { margin-top: 1.2rem; }
  .block-container { max-width: 82rem !important; padding-left: 3rem !important;
    padding-right: 3rem !important; }
</style>
""", unsafe_allow_html=True)


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


# evidence tiers the agent assigns, with a display color + label (Streamlit colored text)
TIER = {
    "ESTABLISHED": ("green", "Established"),
    "EMERGING":    ("orange", "Emerging"),
    "NONE-FOUND":  ("gray", "No published evidence"),
}
_TIER_ORDER = {"ESTABLISHED": 0, "EMERGING": 1, "NONE-FOUND": 2}


def _render_assessments(assessments):
    """Show the agent's cited verdict per candidate, strongest evidence first."""
    for a in sorted(assessments, key=lambda x: _TIER_ORDER.get(x["evidence_tier"], 9)):
        color, label = TIER.get(a["evidence_tier"], ("gray", a["evidence_tier"]))
        st.markdown("**{}** (`{}`) — :{}[{}]".format(
            a["gene"], a["accession"], color, label))
        with st.expander("Agent's cited assessment"):
            st.markdown(a["transcript"])


def _assess_candidate(disease, c):
    """Run the paid cited agent on ONE candidate; record the rate-limit hit and cache the result
    (keyed per disease+accession) so a repeat click never re-spends."""
    st.session_state.hits = st.session_state.get("hits", []) + [time.time()]
    rationale = "physically associates with known {} autoantigens ({})".format(
        DISEASES[disease]["name"], ", ".join(c["connects_to"]))
    prompt = discover.assessment_prompt(c["gene"], c["accession"], DISEASES[disease]["name"], rationale)
    answer, _tools, error = ask_agent(get_agent(disease), prompt)
    if error:
        result = {"evidence_tier": "UNPARSED", "transcript": "The agent hit a problem: " + error}
    else:
        result = {"accession": c["accession"], "gene": c["gene"],
                  "evidence_tier": discover.parse_evidence_tier(answer), "transcript": answer}
    st.session_state.setdefault("assessed", {})[(disease, c["accession"])] = result
    return result


def _tier_badge(verdict):
    color, label = TIER.get(verdict["evidence_tier"], ("gray", verdict["evidence_tier"]))
    st.markdown(":{}[{}]".format(color, label))


def _render_candidate_row(disease, c, precomputed):
    """One live candidate: its evidence, and either a cached/precomputed verdict or an Assess button."""
    verdict = st.session_state.get("assessed", {}).get((disease, c["accession"])) \
        or precomputed.get(c["accession"])
    info, action = st.columns([3, 1])
    with info:
        na = " · nuclear/NA" if c.get("na_associated") else ""
        st.markdown("**{}** (`{}`) — links {} · spec {}{}".format(
            c["gene"], c["accession"], c["partner_support"], c["specificity"], na))
        st.caption("connects to: " + ", ".join(c["connects_to"]))
    with action:
        if verdict:
            _tier_badge(verdict)
        elif st.button("Assess", key="assess_{}_{}".format(disease, c["accession"])):
            if not within_rate_limit():
                st.warning("Hourly limit reached — try later.")
            else:
                with st.spinner("Investigating {}…".format(c["gene"])):
                    verdict = _assess_candidate(disease, c)
                _tier_badge(verdict)
    if verdict and verdict.get("transcript"):
        with st.expander("Agent's cited assessment — {}".format(c["gene"])):
            st.markdown(verdict["transcript"])


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

    # ---- Partner arm: LIVE, steerable epitope-spreading search ----
    st.markdown("#### Partner arm — epitope-spreading neighbors *(live)*")
    if partner:
        with st.expander("How this arm works"):
            st.caption(partner["_metadata"]["method"])

    # STEER: choose which known antigens to spread from (blank = all). Ranking is live and
    # FREE — it runs off a small committed interaction graph, no agent calls.
    pos = discover.positives(disease)
    short_to_acc = {p["short"]: p["accession"] for p in pos}
    picked = st.multiselect("Focus on specific known antigens (leave blank to use all):",
                            list(short_to_acc.keys()), default=[], key="scope_" + disease)
    scope_acc = [short_to_acc[s] for s in picked] or None

    try:
        cands = discover.rank_candidates(disease, scope_accessions=scope_acc)
    except FileNotFoundError:
        st.warning("Live ranking is unavailable (interaction graph not found).")
        cands = []

    # precomputed verdicts (from the saved run) show for free; new ones come from the buttons.
    pre = {a["accession"]: a for a in (partner.get("agent_assessments") if partner else [])}
    if not cands:
        # Explain WHY this scope came up empty instead of a flat "no results". Most often it's
        # because the picked antigen(s) are standalone (e.g. Ro60, a lone RNP with no curated
        # partners), so the epitope-spreading search has nothing to follow — point to antigens
        # that DO form complexes.
        counts = discover.antigen_partner_counts(disease)
        standalone = [s for s in picked if counts.get(s, 0) == 0]
        # prefer a couple of well-known complex antigens that actually have partners; else fall
        # back to whichever known antigens have the most partners.
        suggest = [s for s in ["Sm-B/B'", "La/SSB", "U1-70K", "Ku70"] if counts.get(s, 0) > 0][:2] \
            or [s for s, n in sorted(counts.items(), key=lambda kv: -kv[1]) if n > 0][:2]
        if standalone:
            st.info("**{}** {} no curated interaction partners in UniProt — {} standalone, so the "
                    "epitope-spreading search has nothing to follow. Try a complex-forming antigen "
                    "like **{}** instead.".format(
                        ", ".join(standalone), "has" if len(standalone) == 1 else "have",
                        "it's" if len(standalone) == 1 else "they're", "** or **".join(suggest)))
        else:
            st.caption("No candidates clear the filters for this scope — try adding more antigens.")
    else:
        st.caption("**{}** candidate{} from this scope. *Specificity* = links-to-known-antigens / "
                   "total partners (down-weights promiscuous hubs). Press **Assess** to run the "
                   "cited agent on one (counts against the demo's hourly limit).".format(
                       len(cands), "s" if len(cands) != 1 else ""))
        for c in cands[:12]:
            _render_candidate_row(disease, c, pre)

    # ---- Literature arm: autoantibody-context co-mentions ----
    if literature and literature.get("candidates"):
        st.divider()
        st.markdown("#### Literature arm — autoantibody-context co-mentions")
        with st.expander("How this arm works"):
            st.caption(literature["_metadata"]["method"])
        assessed = {a["accession"]: a for a in literature.get("agent_assessments", [])}
        rows = []
        for c in literature["candidates"]:
            a = assessed.get(c["accession"])
            tier = TIER.get(a["evidence_tier"], ("gray", "—"))[1] if a else "—"
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


# One representative, well-folded structure per disease, shown as a static preview.
HERO = {"lupus": ("P10155", "Ro60"),
        "rheumatoid_arthritis": ("P06733", "Alpha-enolase")}


def _structure_options(disease):
    """[(short, accession)] for autoantigens that have a rendered structure image."""
    opts = []
    for p in DISEASES[disease]["seed"]:
        if p.get("role", "positive") != "positive":
            continue
        if os.path.exists(os.path.join(FIGURES, "{}_alphafold.png".format(p["accession"]))):
            opts.append((p["short"], p["accession"]))
    return opts


def render_structures(disease):
    st.caption("Every structure is an **AlphaFold** prediction, colored by per-residue confidence "
               "(**pLDDT**): blue = very high, orange/red = low (often a flexible / disordered region).")
    # Protein picker → that protein's static render. Static images are reliable; the full WebGL
    # viewer (which misbehaves inside a Streamlit panel) is one click away via the button below.
    opts = _structure_options(disease)
    if opts:
        acc_of = dict(opts)
        hero_acc = HERO.get(disease, (None, None))[0]
        default = next((i for i, (_, a) in enumerate(opts) if a == hero_acc), 0)
        pick = st.selectbox("Pick an autoantigen to view its structure:",
                            [s for s, _ in opts], index=default, key="struct_" + disease)
        acc = acc_of[pick]
        st.image(os.path.join(FIGURES, "{}_alphafold.png".format(acc)), width="stretch",
                 caption="{} ({}) — AlphaFold structure, colored by confidence.".format(pick, acc))
    name = VIEWER.get(disease)
    if not name:
        return
    url = "{}/{}".format(PAGES, name)
    st.link_button("Open the interactive 3D viewer — rotate every structure ↗", url, type="primary")
    st.caption("The full viewer shows all structures together, rotatable and zoomable.")


# ============================================================ sidebar
TOOL_LABEL = {"lookup_protein": "UniProt", "get_structure": "AlphaFold",
              "search_literature": "Europe PMC"}
REPO = "https://github.com/dreamcclure/autoantigen-atlas"
PAGES = "https://dreamcclure.github.io/autoantigen-atlas/figures"

with st.sidebar:
    st.markdown('<p class="side-head">How it works</p>', unsafe_allow_html=True)
    st.markdown(
        "A Claude Opus 4.8 agent with three free tools (UniProt, AlphaFold, Europe PMC) and five "
        "guardrails: it answers **only** from tool results, cites every source, and declines "
        "rather than guesses.")
    st.divider()
    st.markdown('<p class="side-head">Blind-tested</p>', unsafe_allow_html=True)
    st.markdown("Lupus 96% recall · RA 81% · 100% precision & specificity on both.")
    st.divider()
    st.markdown("**Open 3D viewers in a new tab**  \n"
                "[Lupus]({}/structure_viewer.html) · "
                "[RA]({}/rheumatoid_arthritis_viewer.html)".format(PAGES, PAGES))
    st.markdown("[Source on GitHub]({})".format(REPO))
    st.caption("Research aide for scientists — not medical advice.")

# ============================================================ main
# Header: title on the left, the "research aide" blurb on the right.
hc1, hc2 = st.columns([3, 2], vertical_alignment="center")
with hc1:
    st.markdown(
        '<h1 class="atlas-title"><span class="l1">Autoantigen</span>'
        '<span class="l2">Atlas</span></h1>',
        unsafe_allow_html=True,
    )
with hc2:
    st.markdown(
        '<p class="hero-desc">An AI research aide that finds the proteins the immune system '
        'attacks in autoimmune disease — and <b>cites every claim</b>.</p>',
        unsafe_allow_html=True,
    )

# Tab bar, with a SCROLL cue to its right that jumps to the disease section below.
tabs_col, scroll_col = st.columns([8, 1], vertical_alignment="top")
with tabs_col:
    tab_chat, tab_discovery, tab_structures = st.tabs(
        ["Ask the Atlas", "Discovery", "3D Structures"])
with scroll_col:
    st.markdown('<a class="scroll-link" href="#disease-anchor">Scroll ↓</a>',
                unsafe_allow_html=True)

# ---------------------------------------------------------------- tab 1: chat
with tab_chat:
    # Disease selection — moved out of the sidebar into the main flow.
    st.markdown('<div id="disease-anchor"></div>', unsafe_allow_html=True)
    # Disease selection — a black "hero" band with white type (matches dreamcclure.com).
    with st.container(key="disease_band"):
        st.markdown('<p class="band-head">Autoimmune Disease</p>', unsafe_allow_html=True)
        disease = st.selectbox("Autoimmune disease", list(DISEASES.keys()),
                               format_func=lambda k: DISEASES[k]["name"],
                               label_visibility="collapsed")

    # reset the conversation when the disease changes
    if st.session_state.get("disease") != disease:
        st.session_state.disease = disease
        st.session_state.messages = []

    st.caption("Investigating **{}**. Ask whether a protein is an autoantigen, what its structure "
               "looks like, or what the literature says — every answer comes with citations."
               .format(DISEASES[disease]["name"]))

    st.markdown('<p class="sec-head">Ask About a Protein</p>', unsafe_allow_html=True)

    if not st.session_state.get("messages"):
        st.markdown('<p class="tryasking">Suggested Questions</p>', unsafe_allow_html=True)
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

    prompt = st.chat_input("Protein") or st.session_state.pop("pending", None)

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
