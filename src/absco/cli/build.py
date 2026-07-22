"""``absco-build`` — developer Fortran-compile CLI (implemented in a later plan item).

The compilation logic currently lives in :mod:`absco.build_models`; this entry point
will wrap it once the packaging refactor of that module is done (see
docs/packaging_plan.md, item 4).
"""


def main():
    raise SystemExit(
        "absco-build is not implemented yet (see docs/packaging_plan.md, item 4)."
    )


if __name__ == "__main__":
    main()
