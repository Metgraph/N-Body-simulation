#include <stdio.h>
#include <stdlib.h>
#include "stack.h"

stack *new_stack(int size) {
    stack *s;
    s = malloc(sizeof(stack));
    s->array = malloc(size * sizeof(int));
    if (!s->array) {
        printf("Error while creating new stack\n");
        exit(1);
    }
    s->top = -1;
    s->size = size;
    return s;
}

void destroy_stack(stack *s) {
    free(s->array);
    free(s);
}

int is_empty(stack *s) {
    return s->top == -1;
}

int is_full(stack *s) {
    return s->top + 1 == s->size;
}

int size(stack *s) {
    return s->top + 1;
}

void push(stack *s, int v) {
    if (s->top != s->size)
        s->array[++s->top] = v;
}

int pop(stack *s) {
    return s->array[s->top--];
}

int seek(stack *s) {
    return s->array[s->top];
}
