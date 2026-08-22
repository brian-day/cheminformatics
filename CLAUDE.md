# Project context

## Goals

This is a portfolio/demonstration project showcasing cheminformatics and drug-design
concepts (molecule I/O, descriptor calculation, fingerprinting/similarity search). It is
not intended for production use, but code quality, tests, and documentation should still
reflect what a real tool would look like, since it's meant to be shown to others.

The front end can start basic, but I would like to make a medium complexity UI at some point.
For heavier calculations, I'd like the possibility of multi-threaded jobs. 
This project won't be linked to any custom databases yet, but I like to eventually create an API for working
with common public cheminformatics databases.

## Developer background

Developer is a strong software engineer with formal chemistry/physics background, but 
less knowledge of cheminformatics specifically. User wants to use this project also as
a learning opportunity, so err on over-explaining, rather than under-explaining. Remember,
they do have formal training, so more advanced-style explaination are good.
When domain-specific rationale is non-obvious (e.g. why a descriptor like LogP or TPSA
matters, why a particular fingerprint or similarity metric is standard), explain the
"why," not just the "what."

Standard programming concepts don't need extra explanation. User is most experienced in backend
development, but has some knowledge of devops, and a bit less of front end.
Known software languages include python, javascript/typescript, and rust (in decreasing familiarity order).
Has also used matlab and mathematica.
