# STAT 547 Final Project - Persistent Random Walks
This is a collection of scripts used for the exploration of persistent random
walks inspired by their application in Non-Reversible Simulated Tempering
(NRST).  For more background on NRST see Biron-Lattes et al
(https://www.tandfonline.com/doi/full/10.1080/01621459.2024.2335587)
## Getting started
Install the required python libraries. It is advised to do so in a virtual
environment as follows
```
python3 -m venv .env
source .env/bin/activate
pip3 install -r requirements.txt
```
You will need to source the activate script every time you run the project, but
libraries only need to be installed once. The activate script may have different
names on different operating systems. Use the one that is appropriate for your
system.
## Code summary
The following files are available:
- `random_walk.py`: simulate a persistent random walk