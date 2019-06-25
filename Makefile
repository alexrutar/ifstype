.PHONY: test draw profile
test:
	-python3 test.py
draw: test
	-open example.pdf
profile:
	-python3 -m cProfile -o output.pstats test.py
	-gprof2dot -f pstats output.pstats | dot -Tsvg -o output.svg
	-open output.svg
