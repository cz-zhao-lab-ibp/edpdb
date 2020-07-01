#include <stdio.h>
#include <readline/readline.h>
#include <curses.h>

extern int rl_done;

static
int
accept_line(int count, int key)
{
	rl_done = 1;
	return 1;
}

int
smg_create_virtual_keyboard__(int *kbid)
{
	rl_bind_key(RETURN,accept_line);

	return 1;
}

int
smg_read_composed_line__(int *kbid,
                        char *ignored1,
                        char *resultant_string,
                        char *prompt_string,
                        int *prompt_length,
                        int *resultant_length)
{
	char *rline;

	prompt_string[*prompt_length] = '\0';
	
	rline = readline(prompt_string);
	if (rline && *rline)
		add_history(rline);

	strcpy(resultant_string,rline);
	*resultant_length = strlen(rline);

	free(rline);
	printf("\n"); 
	return 1;
}

int
lines_()
{
/*	static int rows, cols;
	rl_get_screen_size(&rows, &cols);
	printf("%d %d \n",LINES,COLS ); */
	return LINES;
}
