.PHONY: clean

app: app.c project.o result.h panic.h

test_project: test_project project.o

project.o: project.c project.h

clean:
	rm -f app test_project project.o
