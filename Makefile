.PHONY: test draw profile
test:
	-python3 test.py
draw: test
	-open example.pdf
profile:
	-python3 -m cProfile -o profile/output.pstats test.py
	-gprof2dot -f pstats profile/output.pstats | dot -Tsvg -o profile/output.svg
	-open profile/output.svg
