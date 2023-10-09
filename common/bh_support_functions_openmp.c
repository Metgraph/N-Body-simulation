void s_center_of_mass(Octtree *tree, Octnode *node) {
    Octnode *child;
    double new_mass;

    if (node->ents == 1)
        return;
    for (int n = 0; n < 8; n++) {
        if (node->children[n] != -1)
            s_center_of_mass(tree, &tree->nodes[node->children[n]]);
    }
    for (int n = 0; n < 8; n++) {
        if (node->children[n] != -1) {
            child = &tree->nodes[node->children[n]];

            new_mass = node->mass + child->mass;

            node->center.x = (child->center.x * child->mass / new_mass) + (node->center.x * node->mass / new_mass);
            node->center.y = (child->center.y * child->mass / new_mass) + (node->center.y * node->mass / new_mass);
            node->center.z = (child->center.z * child->mass / new_mass) + (node->center.z * node->mass / new_mass);

            node->mass = new_mass;
        }
    }
}

void s_get_bounding_box(Entity ents[], int ents_sz, double *max) {
    *max = 0.0;

    for (int i = 0; i < ents_sz; i++) {
        *max = fabs(ents[i].pos.x) > *max ? fabs(ents[i].pos.x) : *max;
        *max = fabs(ents[i].pos.y) > *max ? fabs(ents[i].pos.y) : *max;
        *max = fabs(ents[i].pos.z) > *max ? fabs(ents[i].pos.z) : *max;
    }
    //printf("Max: %lf\n", *max);
    *max *= 2;
}

void s_get_energy(Entity *ents, int sz, double *KE, double *PE) {

    // calculate KE
    double local_KE=0.0;
    for (int i = 0; i < sz; i++) {
        RVec3 vel = ents[i].vel;
        double mass = ents[i].mass;

        local_KE += vel.x * vel.x *mass;
        local_KE += vel.x * vel.x *mass;
        local_KE += vel.x * vel.x *mass;
    }
    *KE=local_KE*0.5;

    //calculate PE
    double local_PE = 0.0;

    //calculate only the upper triangle in matrix instead of calculate whole matrix and take only the
    //upper triangle values like using np.triu
    for (int i = 0; i < sz; i++) {
        RVec3 *i_pos = &ents[i].pos;
        for (int j = i; j < sz; j++) {
            double dx, dy, dz, D;
            RVec3 *j_pos = &ents[j].pos;

            dx = i_pos->x - j_pos->x;
            dy = i_pos->y - j_pos->y;
            dz = i_pos->z - j_pos->z;
            D = sqrt(dx * dx + dy * dy + dz * dz);
            D = D > 0 ? 1.0 / D : D;

            local_PE += -(ents[i].mass * ents[j].mass) * D;
        }
    }
    *PE=local_PE*BIG_G;
}


