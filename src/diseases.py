"""
The platform's disease registry.
================================
This is the keystone of the "platform" idea: a disease is just *data* here, not code.
To support a new autoimmune disease, you add an entry to DISEASES — the fetch, build,
and evaluate steps all read from this one place. Lupus was the flagship; rheumatoid
arthritis proves the same engine generalizes.

Each seed entry: short name, UniProt accession, gene symbol, the clinical antibody
test, and a role ("positive" = established autoantigen, "negative_control" = an
abundant protein NOT targeted in that disease, used to test specificity).

Every accession was verified live against UniProt (lupus 2026-11-01, RA 2026-11-02).
"""

DISEASES = {
    "lupus": {
        "name": "Systemic lupus erythematosus (SLE)",
        "source": "Autoantibodies in SLE: Diagnostic and Pathogenic Insights (2025), PMC12386345",
        # Disease-specific reasoning the agent uses to judge autoantigen status.
        # Keeping it here (data) rather than in the agent code keeps the platform generic.
        "guidance": (
            "In SLE, a protein that is an established core component of a recognized autoantigen "
            "complex — the Sm core / snRNP, U1 snRNP, nucleosome (histones), ribosomal P, Ku, or "
            "the Ro/La particle — and is captured by that complex's named clinical antibody test "
            "(anti-Sm, anti-U1-RNP, anti-histone/nucleosome, anti-ribosomal P, anti-Ku, anti-Ro/La) "
            "is itself an established SLE autoantigen. Complex-level evidence plus the documented "
            "membership is sufficient; you need not find protein-chain-specific papers. This does NOT "
            "apply to proteins that merely interact with, or are unrelated to, those complexes."
        ),
        "seed": [
            # --- Ro / La ---
            {"short": "Ro60", "accession": "P10155", "gene": "RO60", "antibody_test": "anti-Ro/SSA", "role": "positive"},
            {"short": "Ro52", "accession": "P19474", "gene": "TRIM21", "antibody_test": "anti-Ro52", "role": "positive"},
            {"short": "La/SSB", "accession": "P05455", "gene": "SSB", "antibody_test": "anti-La/SSB", "role": "positive"},
            # --- Sm ---
            {"short": "Sm-B/B'", "accession": "P14678", "gene": "SNRPB", "antibody_test": "anti-Sm", "role": "positive"},
            {"short": "Sm-D1", "accession": "P62314", "gene": "SNRPD1", "antibody_test": "anti-Sm", "role": "positive"},
            {"short": "Sm-D2", "accession": "P62316", "gene": "SNRPD2", "antibody_test": "anti-Sm", "role": "positive"},
            {"short": "Sm-D3", "accession": "P62318", "gene": "SNRPD3", "antibody_test": "anti-Sm", "role": "positive"},
            # --- U1-RNP ---
            {"short": "U1-70K", "accession": "P08621", "gene": "SNRNP70", "antibody_test": "anti-U1-RNP", "role": "positive"},
            {"short": "U1-A", "accession": "P09012", "gene": "SNRPA", "antibody_test": "anti-U1-RNP", "role": "positive"},
            {"short": "U1-C", "accession": "P09234", "gene": "SNRPC", "antibody_test": "anti-U1-RNP", "role": "positive"},
            # --- nucleosome / histones ---
            {"short": "Histone H4", "accession": "P62805", "gene": "H4C1", "antibody_test": "anti-histone/nucleosome", "role": "positive"},
            {"short": "Histone H3.1", "accession": "P68431", "gene": "H3C1", "antibody_test": "anti-histone/nucleosome", "role": "positive"},
            {"short": "Histone H2B", "accession": "P62807", "gene": "H2BC4", "antibody_test": "anti-histone/nucleosome", "role": "positive"},
            {"short": "Histone H2A", "accession": "P0C0S8", "gene": "H2AC11", "antibody_test": "anti-histone/nucleosome", "role": "positive"},
            # --- ribosomal P ---
            {"short": "Ribosomal P0", "accession": "P05388", "gene": "RPLP0", "antibody_test": "anti-ribosomal P", "role": "positive"},
            {"short": "Ribosomal P1", "accession": "P05386", "gene": "RPLP1", "antibody_test": "anti-ribosomal P", "role": "positive"},
            {"short": "Ribosomal P2", "accession": "P05387", "gene": "RPLP2", "antibody_test": "anti-ribosomal P", "role": "positive"},
            # --- complement C1q ---
            {"short": "C1q A", "accession": "P02745", "gene": "C1QA", "antibody_test": "anti-C1q", "role": "positive"},
            {"short": "C1q B", "accession": "P02746", "gene": "C1QB", "antibody_test": "anti-C1q", "role": "positive"},
            {"short": "C1q C", "accession": "P02747", "gene": "C1QC", "antibody_test": "anti-C1q", "role": "positive"},
            # --- DNA-repair (overlap) ---
            {"short": "Ku70", "accession": "P12956", "gene": "XRCC6", "antibody_test": "anti-Ku", "role": "positive"},
            {"short": "Ku80", "accession": "P13010", "gene": "XRCC5", "antibody_test": "anti-Ku", "role": "positive"},
            # --- proliferation ---
            {"short": "PCNA", "accession": "P12004", "gene": "PCNA", "antibody_test": "anti-PCNA", "role": "positive"},
            # --- antiphospholipid ---
            {"short": "Beta-2-GPI", "accession": "P02749", "gene": "APOH", "antibody_test": "anti-beta2GPI", "role": "positive"},
            # --- negative controls (NOT established lupus autoantigens) ---
            {"short": "GAPDH", "accession": "P04406", "gene": "GAPDH", "antibody_test": "(none)", "role": "negative_control"},
            {"short": "Beta-actin", "accession": "P60709", "gene": "ACTB", "antibody_test": "(none)", "role": "negative_control"},
            {"short": "Serum albumin", "accession": "P02768", "gene": "ALB", "antibody_test": "(none)", "role": "negative_control"},
            {"short": "Hemoglobin beta", "accession": "P68871", "gene": "HBB", "antibody_test": "(none)", "role": "negative_control"},
            {"short": "Tubulin beta", "accession": "P07437", "gene": "TUBB", "antibody_test": "(none)", "role": "negative_control"},
        ],
    },
    "rheumatoid_arthritis": {
        "name": "Rheumatoid arthritis (RA)",
        "source": "Autoantibodies in Rheumatoid Arthritis (2024), PMC11139241",
        # NOTE: most RA antigens are the *citrullinated* (post-translationally modified)
        # form; UniProt/AlphaFold hold the unmodified protein. Documented, like dsDNA for lupus.
        "guidance": (
            "In RA, the two main antibody systems are ACPA (anti-citrullinated protein antibodies) "
            "and rheumatoid factor (RF). ACPA recognize proteins whose arginine residues have been "
            "citrullinated (a post-translational modification); RF targets the Fc (constant) region "
            "of immunoglobulin G. Judge a protein an established RA autoantigen only if your tools "
            "provide evidence that its (citrullinated) form is an established ACPA target, that it "
            "is the RF target, or that it is otherwise an established autoantibody target in RA. "
            "Because UniProt and AlphaFold hold the unmodified protein, do not reject a protein "
            "merely because the database entry is unmodified — but do require tool evidence linking "
            "it to RA autoimmunity before answering YES. Additionally, a protein that your tools "
            "show is a chain, subunit, or member of a larger protein or protein family that your "
            "tools establish as an RA autoantibody target is itself an established RA autoantigen "
            "(e.g. an individual chain of a multi-chain target protein, or an individual member of a "
            "protein family whose citrullinated form is targeted) — provided the tools establish "
            "that the parent protein or family is the RA target; do not assume it."
        ),
        "seed": [
            # --- ACPA (anti-citrullinated protein) targets ---
            {"short": "Vimentin", "accession": "P08670", "gene": "VIM", "antibody_test": "ACPA (anti-MCV)", "role": "positive"},
            {"short": "Alpha-enolase", "accession": "P06733", "gene": "ENO1", "antibody_test": "ACPA (anti-CEP-1)", "role": "positive"},
            {"short": "Fibrinogen alpha", "accession": "P02671", "gene": "FGA", "antibody_test": "ACPA", "role": "positive"},
            {"short": "Fibrinogen beta", "accession": "P02675", "gene": "FGB", "antibody_test": "ACPA", "role": "positive"},
            {"short": "Fibrinogen gamma", "accession": "P02679", "gene": "FGG", "antibody_test": "ACPA", "role": "positive"},
            {"short": "Fibronectin", "accession": "P02751", "gene": "FN1", "antibody_test": "ACPA", "role": "positive"},
            {"short": "Type II collagen", "accession": "P02458", "gene": "COL2A1", "antibody_test": "ACPA / anti-CII", "role": "positive"},
            {"short": "Filaggrin", "accession": "P20930", "gene": "FLG", "antibody_test": "ACPA (anti-filaggrin)", "role": "positive"},
            {"short": "Histone H4", "accession": "P62805", "gene": "H4C1", "antibody_test": "ACPA (cit-histone)", "role": "positive"},
            {"short": "Histone H2B", "accession": "P62807", "gene": "H2BC4", "antibody_test": "ACPA (cit-histone)", "role": "positive"},
            {"short": "Histone H3.1", "accession": "P68431", "gene": "H3C1", "antibody_test": "ACPA (cit-histone)", "role": "positive"},
            {"short": "Histone H2A", "accession": "P0C0S8", "gene": "H2AC11", "antibody_test": "ACPA (cit-histone)", "role": "positive"},
            # --- other established RA autoantigens ---
            {"short": "BiP/GRP78", "accession": "P11021", "gene": "HSPA5", "antibody_test": "anti-BiP", "role": "positive"},
            {"short": "RA33", "accession": "P22626", "gene": "HNRNPA2B1", "antibody_test": "anti-RA33", "role": "positive"},
            {"short": "Pentraxin 3", "accession": "P26022", "gene": "PTX3", "antibody_test": "anti-PTX3", "role": "positive"},
            {"short": "IgG1 Fc", "accession": "P01857", "gene": "IGHG1", "antibody_test": "rheumatoid factor", "role": "positive"},
            # --- negative controls (NOT established RA autoantigens; beta-actin EXCLUDED — it is one) ---
            {"short": "Hemoglobin beta", "accession": "P68871", "gene": "HBB", "antibody_test": "(none)", "role": "negative_control"},
            {"short": "Insulin", "accession": "P01308", "gene": "INS", "antibody_test": "(none)", "role": "negative_control"},
            {"short": "Lysozyme C", "accession": "P61626", "gene": "LYZ", "antibody_test": "(none)", "role": "negative_control"},
            {"short": "Myoglobin", "accession": "P02144", "gene": "MB", "antibody_test": "(none)", "role": "negative_control"},
            {"short": "Tubulin beta", "accession": "P07437", "gene": "TUBB", "antibody_test": "(none)", "role": "negative_control"},
        ],
    },
}


def seed_for(disease):
    """Return the seed list for a disease key (e.g. 'lupus'), with a clear error if unknown."""
    if disease not in DISEASES:
        raise ValueError("Unknown disease '{}'. Known: {}".format(
            disease, ", ".join(DISEASES)))
    return DISEASES[disease]["seed"]
