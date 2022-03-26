#include "common.h"
#include "assemble.h"
#include "scaffold.h"
#include "gap_close.h"
#include "bubble_map.h"
#include "correct.h"
#include <omp.h>

void print_usage(void)
{
    fputs("Usage: platanus Command [options]\n\nCommand: assemble, scaffold, gap_close, bubble_map correct\n", stderr);
}

void print_version(void)
{
    fprintf(stderr,"Version: %s\n", VERSION);
}

int main(int argc, char **argv)
{
	print_version();
    put_command_log(argc, argv, stderr);

    if (argc > 1) {
        if (strcmp(argv[1], "assemble") == 0) {
            option_assemble_t opt_ass = {};

            option_assemble_init(&opt_ass);
            if (option_assemble_parse(&opt_ass, argc, argv) != 0) {
                print_assemble_usage();
                return 1;
            }
            if (opt_ass.t == 1)
                assemble(&opt_ass);
            else
                assemble_mt(&opt_ass);
            option_assemble_destroy(&opt_ass);
        }
        else if (strcmp(argv[1], "scaffold") == 0) {
            option_scaffold_t opt_sca = {};

            option_scaffold_init(&opt_sca);
            if (option_scaffold_parse(&opt_sca, argc, argv) != 0) {
                print_scaffold_usage();
                return 1;
            }
            scaffold_mt(&opt_sca);
            option_scaffold_destroy(&opt_sca);
        }
        else if (strcmp(argv[1], "gap_close") == 0) {
            option_gap_close_t opt_gap = {};

            option_gap_close_init(&opt_gap);
            if (option_gap_close_parse(&opt_gap, argc, argv) != 0) {
                print_gap_close_usage();
                return 1;
            }
            gap_close_mt(&opt_gap);
            option_gap_close_destroy(&opt_gap);
        }
        else if (strcmp(argv[1], "bubble_map") == 0) {
            option_bubble_map_t opt_bub = {};

            option_bubble_map_init(&opt_bub);
            if (option_bubble_map_parse(&opt_bub, argc, argv) != 0) {
                print_bubble_map_usage();
                return 1;
            }
            bubble_map(&opt_bub);
            option_bubble_map_destroy(&opt_bub);
        }
        else if (strcmp(argv[1], "correct") == 0) {
            option_correct_t opt_cor = {};

            option_correct_init(&opt_cor);
            if (option_correct_parse(&opt_cor, argc, argv) != 0) {
                print_correct_usage();
                return 1;
            }
            if (opt_cor.t == 1)
                correct(&opt_cor);
            else
                correct_mt(&opt_cor);
            option_correct_destroy(&opt_cor);
        }
        else if (!strcmp(argv[1], "-v") || !strcmp(argv[1], "-version") || !strcmp(argv[1], "--version")) {
            print_version();
			return 0;
		}
        else {
            print_usage();
            return 1;
        }
    }
    else {
        print_usage();
        return 1;
    }
    return 0;
}
