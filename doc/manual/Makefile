NAME = manual
FILE = manual.tex
RM = /bin/rm -f
all:	$(NAME) clean
$(NAME): $(FILE)
	-$(RM) $(NAME).pdf
	rubber -d $^
clean:
	rubber --clean $(FILE)
	-$(RM) *~
