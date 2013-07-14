#!/usr/bin/env perl
#
# IsolateTargetTree.pl
# Usage: IsolateTargetTree.pl <target_tree_id> <trees_file>

use strict;
use warnings;

my $target_tree = shift;

# Print up to the first tree declaration, then print only the target tree
my $is_preamble = 1;
while (<>) {
    $is_preamble = 0 if /^tree STATE_.*/;
    print if $is_preamble or /$target_tree\s/;
}

print "End;\n";

