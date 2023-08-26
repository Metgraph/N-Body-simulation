void print_tree_rec(Octtree *tree, int id, char *space, uint depth) {
    Octnode *node = &tree->nodes[id];
    // how much divide
    uint temp = 1 << depth;
    double border = (tree->max) / (double)temp;
    printf("%sid: %d, (x:%lf, y:%lf, z:%lf), border: %lf\n", space, id,
           node->center.x, node->center.y, node->center.z, border);
    if (node->ents > 1) {

        int i;
        for (i = depth * 4; i < depth * 4 + 4; i++) {
            space[i] = ' ';
        }
        space[i] = '\0';
        for (int i = 0; i < 8; i++) {
            if (node->children[i] > -1) {
                print_tree_rec(tree, node->children[i], space, depth + 1);
            }
        }
        space[depth * 4] = '\0';
    }
}

// used for debug
void print_tree(Octtree *tree) {
    uint sz_space = 4 * 40;
    char *space = malloc(sz_space * sizeof(char));
    space[0] = '\0';
    print_tree_rec(tree, tree->root, space, 0);
    free(space);
}

void print_tree_indented(Octtree *tree, Octnode *node, int tabs, int pos) {
    printf("Body position %d -> ", pos);
    if (node->ents == 1) {
        printf("Total ents: %d, parent: %d\n", node->ents, node->parent);
        return;
    }
    printf("Total ents: %d, parent: %d\n", node->ents, node->parent);
    for (int i = 0; i < 8; i++) {
        if (node->children[i] == -1) {
            for (int j = 0; j < tabs; j++)
                printf("\t|");
            printf("Child %d is empty\n", i);
        } else {
            for (int j = 0; j < tabs; j++)
                printf("\t|");
            printf("Child %d: ", i);
            print_tree_indented(tree, &tree->nodes[node->children[i]], tabs + 1,
                                node->children[i]);
        }
    }
}

