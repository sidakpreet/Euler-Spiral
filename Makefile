CFLAGS = -g -DHAVE_CONFIG_H -fPIC
CC = g++

# Colour coding of the output text
NO_COLOR	  = \033[0m
WHITE_COLOR   = \033[37;01m
OK_COLOR      = \033[32;01m
ERROR_COLOR   = \033[31;01m
WARN_COLOR    = \033[33;01m
MAGENTA_COLOR = \033[35;01m


# Make all

all:	test
	@echo "$(MAGENTA_COLOR)*** ALL Files successfully created!$(NO_COLOR)"


test: test.cpp
	@echo "$(WHITE_COLOR)* Creating Test file$(NO_COLOR)"
	@$(CC) $(INCLUDES) $(CFLAGS) $^ -o $@ $(LIBS) $(opencv_LIBS)
	@echo "$(WHITE_COLOR)*$(OK_COLOR) .......................................[OK]$(NO_COLOR)"
	
clean:
	@rm test
	@echo "$(WARN_COLOR)* All Files cleared!!$(NO_COLOR)"
