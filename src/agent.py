"""
Phase 2 — The Agent
===================
A LangGraph agent, powered by Claude, that researches lupus autoantigens using our
three tools (UniProt protein lookup, AlphaFold structure, Europe PMC literature) and
writes a CITED summary for a human scientist to review.

Run it from the project folder:
    python3 src/agent.py "What is Ro60 and why is it a lupus autoantigen?"
Or with no question, for an interactive chat:
    python3 src/agent.py

Needs ANTHROPIC_API_KEY in a .env file at the project root (see .env.example).
"""

import os
import sys

from dotenv import load_dotenv

# --- Load the secret API key from .env at the project root -------------------
HERE = os.path.dirname(os.path.abspath(__file__))              # .../src
PROJECT_ROOT = os.path.dirname(HERE)                            # project root
load_dotenv(dotenv_path=os.path.join(PROJECT_ROOT, ".env"))

from langchain_anthropic import ChatAnthropic
from langgraph.prebuilt import create_react_agent

from tools import ALL_TOOLS

# The model — Anthropic's most capable, best at multi-step tool-using research.
MODEL = "claude-opus-4-8"

# --- The guardrails, written as the agent's standing instructions ------------
SYSTEM_PROMPT = """\
You are a research assistant that helps scientists study autoantigens in lupus
(systemic lupus erythematosus, SLE) — the self-proteins the immune system mistakenly
attacks. You have three tools: lookup_protein (UniProt facts), get_structure
(AlphaFold 3D-structure confidence), and search_literature (Europe PMC papers).

Follow these rules without exception:

1. CITE OR DON'T CLAIM. State a fact only if a tool returned it in this conversation,
   and attach the source the tool gave you (a UniProt, AlphaFold, or paper URL) right
   next to the claim. If you have no tool evidence for something, say so plainly.

2. NO HALLUCINATION. Do not rely on your own prior knowledge of biology for factual
   claims — use the tools. Your memory may be wrong or outdated; the databases are not.

3. HUMAN IN THE LOOP. You are an aide to an expert, not an autonomous decision-maker.
   Present findings for the scientist to review and do not overstate certainty. When
   evidence is thin or mixed, say so.

4. STAY IN SCOPE. Help only with lupus autoantigens and the surrounding protein biology.
   Politely decline unrelated requests.

5. FAIL HONESTLY. If a tool returns an error or no data, say that clearly and suggest a
   next step — never invent an answer to fill the gap.

Typical workflow for a protein question: call lookup_protein for the facts, then
get_structure (using the accession it returns) for structural confidence, then
search_literature for published evidence. Finish with a short, clearly-cited summary.

Judging autoantigen status (based only on tool evidence):
- A protein that lookup_protein shows is an established component of a recognized SLE
  autoantigen complex — the Sm core / snRNP, U1 snRNP, nucleosome (histones), ribosomal
  P, Ku, or Ro/La particle — and is captured by that complex's named clinical antibody
  test (anti-Sm, anti-U1-RNP, anti-histone/nucleosome, anti-ribosomal P, anti-Ku,
  anti-Ro/La) is itself an established SLE autoantigen. Complex-level evidence plus the
  documented membership is sufficient; you do not need protein-chain-specific papers.
  (This does NOT apply to proteins that merely interact with, or are unrelated to, those
  complexes — membership must be as a recognized core component of the autoantigen.)
- If a literature search returns off-topic results, refine the query with explicit
  disease and antibody terms (e.g., "anti-PCNA antibody systemic lupus erythematosus")
  before concluding there is no evidence.
"""


def build_agent():
    """Wire Claude + the guardrails + the three tools into a ready-to-use agent."""
    model = ChatAnthropic(
        model=MODEL,
        max_tokens=8192,                   # room for internal thinking + the answer
        thinking={"type": "adaptive"},     # Claude decides how hard to think, per step
    )
    return create_react_agent(model, ALL_TOOLS, prompt=SYSTEM_PROMPT)


def _text(content):
    """Pull the plain-text answer out, whether it's a string or a list of blocks."""
    if isinstance(content, list):
        return "".join(
            block.get("text", "")
            for block in content
            if isinstance(block, dict) and block.get("type") == "text"
        )
    return content


def ask(agent, question):
    """Send one question to the agent and return its final cited answer."""
    result = agent.invoke({"messages": [{"role": "user", "content": question}]})
    return _text(result["messages"][-1].content)


def main():
    agent = build_agent()

    # One-shot mode: a question was passed on the command line.
    if len(sys.argv) > 1:
        question = " ".join(sys.argv[1:])
        print("\nQ: {}\n".format(question))
        print(ask(agent, question))
        return

    # Interactive mode: chat until the user types quit.
    print("Lupus Autoantigen Assistant (Claude Opus 4.8). Type 'quit' to exit.\n")
    while True:
        try:
            question = input("You: ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break
        if question.lower() in ("quit", "exit", "q", ""):
            break
        print("\nAssistant:", ask(agent, question), "\n")


if __name__ == "__main__":
    main()
