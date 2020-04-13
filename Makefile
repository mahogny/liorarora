putim:
	scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" analyze_images_b* mahogny@ebi-login.ebi.ac.uk:/homes/mahogny/common/liora/microscopy/rora

getim:
	scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" mahogny@ebi-login.ebi.ac.uk:/homes/mahogny/common/liora/microscopy/rora/images.out/* images.out/

putbed:
	scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" *.bed mahogny@ebi-login.ebi.ac.uk:/homes/mahogny/common/liora/

getfimo:
	scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" mahogny@ebi-login.ebi.ac.uk:/homes/mahogny/common/liora/fimo_upstream/fimo.txt   out/upstream.csv
	scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" mahogny@ebi-login.ebi.ac.uk:/homes/mahogny/common/liora/fimo_downstream/fimo.txt out/downstream.csv


summarize:
	cat images.out/*.meta* > images.meta.csv

foo:
	scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" ../upload_cyto.sh  mahogny@ebi-login.ebi.ac.uk:/homes/mahogny/common/data/henriksson_20180530_cytoscreen/
