.PHONY: run test clean lint style

run:
	Rscript DGE_pipeline.R

test:
	Rscript test/run_tests.R

lint:
	Rscript -e "lintr::lint_dir('.')"

style:
	Rscript -e "styler::style_dir('.')"

clean:
	rm -rf Output/*
