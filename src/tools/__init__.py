"""The agent's toolbox. Each tool is a function the AI is allowed to call.

ALL_TOOLS is the list we hand to the agent in src/agent.py.
"""

from .protein import lookup_protein
from .structure import get_structure
from .literature import search_literature

ALL_TOOLS = [lookup_protein, get_structure, search_literature]
