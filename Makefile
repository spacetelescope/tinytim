SRC = src
.PHONY: all
all:
	$(MAKE) -C $(SRC) all

install:
	$(MAKE) -C $(SRC) install

clean:
	$(MAKE) -C $(SRC) clean

