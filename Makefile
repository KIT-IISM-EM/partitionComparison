

all: install

clean:
	rm -rf ./partitionComparison.Rcheck
	rm -f partitionComparison_*.tar.gz

build: clean
	R CMD build partitionComparison

install: clean
	R CMD INSTALL partitionComparison

check: build
	R CMD check partitionComparison_*.tar.gz

cran: build
	R CMD check --as-cran partitionComparison_*.tar.gz

doc: 
	rm -f partitionComparison.pdf
	R -e 'library("devtools"); devtools::document("partitionComparison")'
	R CMD Rd2pdf partitionComparison
