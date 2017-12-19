

all: install

clean:
	rm -rf ./partitionComparison.Rcheck
	rm -f partitionComparison_*.tar.gz
	rm -f partitionComparison.pdf

build: clean
	R CMD build partitionComparison

install: clean
	R CMD INSTALL partitionComparison

check: clean
	R cmd check partitionComparison --as_cran

doc: clean
	R CMD Rd2pdf partitionComparison
