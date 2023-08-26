#ifndef STACK
#define STACK

typedef struct S {
    int top;
    int size;
    int *array;
} stack;

stack *new_stack(int size);
void destroy_stack(stack *s);
int is_empty(stack *s);
int is_full(stack *s);
int size(stack *s);
void push(stack *s, int v);
int pop(stack *s);
int seek(stack *s);

#endif
