.PHONY: test draw
test:
	-python3 -m ifstype
draw: test
	-open example.pdf
