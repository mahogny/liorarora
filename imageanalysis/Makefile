get:
	scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" mahogny@ebi-login.ebi.ac.uk:/homes/mahogny/common/liora/microscopy/rora/images.out/* images.out/

collect:
	cat images.out/*meta* > meta.out.csv
