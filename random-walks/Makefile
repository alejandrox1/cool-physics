IMG := walk
WORKDIR := /app

CXX = g++
CXXFLAGS = -Wall -Wextra -Werror
LIBS = -lgsl
EXEC = main.exec
OBJECTS = main.o

.PHONY: all
all: dev

.PHONY: compile
compile: $(OBJECTS)
	$(CXX) -o $(EXEC) $^ $(LIBS)

%.o: %.cpp
	$(CXX) -o $@ -c $^ $(CXXFLAGS)

.PHONY: clean
clean:
	rm -rf $(OBJECTS)

.PHONY: build
build:
	docker build -t $(IMG) .

.PHONY: dev
dev: build
	docker run --rm -it -v $(PWD):$(WORKDIR) $(IMG) bash

.PHONY: run
run: build
	docker run --rm -it -v $(PWD):$(WORKDIR) $(IMG)
