"""PRISM: CLI Entry Point.

Delegates to the CLI controller so the pipeline can be invoked as a module via
``python -m PRISM``.
"""

from PRISM.controller import main

if __name__ == "__main__":
    main()
