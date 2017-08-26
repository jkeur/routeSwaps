/** Algorithm simpleSort to get an optimal cycle decomposition of the move graph G'
 * belonging to the problem of routing qubits from distribution x to y
 * in the fully connected star graph G = (V, E).
 *
 * Created by J. Keur
 * August 2017
 */

//#define PRINT_SWAPS     // Print swaps done
//#define PRINT_NUM       // Print qubit no.

extern void setW();
extern char alg;            // Solver algorithm: ROUTE_SIMPLE or ROUTE_SWAPS

/** Test if node is centre node
 * 170410 Created
 */
unsigned isC(const unsigned j)
{
    return ((j % (m + 1)) == 0) ? 1 : 0;
}

/** Take action: swap over edge ei
 * Returns non-zero value if succeeded, returns 0 if failed
 * 170411 Created
 */
char swap(const unsigned i, const unsigned j)
{
    const unsigned tmp = x[i];
    const unsigned gi = i/(m+1);
    const unsigned gj = j/(m+1);

    // If the swap cannot be done now
    if ( ((c2use[gi] & IGNORE_C) && (c2use[gi] != CORRECT))
            || (isC(j) && ((c2use[gj] & IGNORE_C) && (c2use[gj] != CORRECT))) )
        return 0;

    // Swap qubits
    x[i] = x[j];
    x[j] = tmp;
#ifdef PRINT_SWAPS
    // Print the swap done
    if (isC(j))
    {
        printf("%*u", Nd, i+1);
#ifdef PRINT_NUM
        COLOR_NUM;
        printf("(%*u)", Nd, x[i]+1);
#endif // PRINT_NUM
        COLOR_TEXT;         // 170417
        printf("-");
        NORMAL_TEXT;        // Standard color
#ifdef PRINT_NUM
        printf("%*u", Nd, j+1);
        COLOR_NUM;
        printf("(%*u) ", Nd, x[j]+1);
        NORMAL_TEXT;
#else
        printf("%*u ", Nd, j+1);
#endif // PRINT_NUM
    }
    else
    {
#ifdef PRINT_NUM
        printf("%*u", Nd, i+1);
        COLOR_NUM;
        printf("(%*u)", Nd, x[i]+1);
        NORMAL_TEXT;
        printf("-%*u", Nd, j+1);
        COLOR_NUM;
        printf("(%*u) ", Nd, x[j]+1);
        NORMAL_TEXT;
#else
        printf("%*u-%*u ", Nd, i+1, Nd, j+1);
#endif
    }
#endif // PRINT_SWAPS
    if (isC(j))
    {
        Nsb++;                      // One expensive swap done
    }

    c2use[i/(m+1)] = BEING_USED;    // This centre is swapped now; it cannot be used in this stage
    c2use[j/(m+1)] = BEING_USED;    // This centre is swapped now
    Ns++;                           // Increase swap counter
    setW();
    if (isC(j))                     // If a swap between centres was done
    {
        if (W[i/(m+1)][i/(m+1)] == m + 1) // If all numbers are in group Gi
        {
            if (isC(x[i]))          // If this centre has the right number
                c2use[i/(m+1)] = SORTED;
            else
                c2use[i/(m+1)] = CORRECT | BEING_USED;
        }
        if (W[j/(m+1)][j/(m+1)] == m + 1) // If all numbers are in group Gj
        {
            if (isC(x[j]))          // If this centre has the right number
                c2use[j/(m+1)] = SORTED;
            else
                c2use[j/(m+1)] = CORRECT | BEING_USED;
        }
    }

    return OK;      // Swap done
}

/** If possible, do a swap with Gi such that Gi gets a qubit from Gj
 * 170520 Created
 * /
static void doS(const unsigned gi, const unsigned gj)
{
    if ( (D[gi] != gi) && (D[gj] == gi) )       // If a swap is possible
        swap(gi*(m+1), gj*(m+1));               // Do the swap
}*/

/** Set xi such that it should be moved to group Gj.
 * Return 1 if a swap is done
 * Return 0 if the centre already has a number for group Gj
 * 170520 Created
 */
static char setN(const unsigned gi, const unsigned gj)
{
    unsigned l, di;
    const unsigned i = gi*(m+1);

    di = getDestStar(i);            // Get destination star
    if (di == gj)
        return 0;                   // If centre already has a number for Gj
    for (l = 1; l <= m; l++)        // For each leaf
        if (getDestStar(i + l) == gj) // If x(i+l) should be moved to Gj
        {
            swap(i, i+l);           // Swap the leaf number
            return 1;               // Return nonzero; a swap is done
        }
    return -2;  // ERROR
}

/** Ensure that the number xi in the centre of star gi should be moved outwards.
 * Return 1 if a swap is done
 * Return 0 if the centre already has a number for group Gj != Gi
 * 170520 Created
 */
static char setOut(const unsigned gi)
{
    unsigned l, di;
    const unsigned i = gi*(m+1);

    di = getDestStar(i);            // Get destination star
    if (di != gi)
        return 0;                   // If centre already has a number for Gj != Gi
    for (l = 1; l <= m; l++)        // For each leaf
        if (getDestStar(i + l) != gi) // If x(i+l) should be moved outwards
        {
            swap(i, i+l);           // Swap the leaf number
            return 1;               // Return nonzero; a swap is done
        }
    return -2;  // ERROR
}

/** This is a trivial sorting algorithm
 * which firstly places all numbers in the first group, then in the second, etc.
 * 170520 Created
 */
void routeSimple()
{
    unsigned gi, gj, si, sj;

#ifdef PRINT_SWAPS
    puts("RouteSimple is solving the problem...");
#endif // PRINT_SWAPS

    newRound();

    for (gi = 0; gi < k; gi++)      // For each group Gi
    {
        //printf("i_%u ", gi+1);
        gj = gi;
        // Get the right numbers into group Gi
        for (gj = gi+1; gj < k; gj++) // For each other group Gj > Gi
        {
            if ( (gj == gi) || (W[gj][gi] == 0) ) // If Gj!=Gi does not have any number for group Gi
                continue;
            //printf("j_%u ", gj+1);
            newStage();
            si = setOut(gi);                // si = 1 if a swap in Gi is done
            sj = setN(gj, gi);              // sj = 1 if a swap in Gj is done
            if ((si == 1) || (sj == 1))     // If >=1 swap is done such that D(xi)!=Gi and D(xj) = Gi
                newStage();                 // Begin a new stage
            //doS(gi, gj);            // Do the swap
            swap(gi*(m+1), gj*(m+1));       // Do the swap
            while (W[gj][gi])               // While Gj has numbers for Gi
            {
                newStage();                 // Bring a number from Gj to Gi
                setOut(gi);                 // Swap such that dest(xi) != Gi
                setN(gj, gi);               // Swap such that dest(xj) = Gi
                newStage();                 // Begin a new stage
                //doS(gi, gj);
                swap(gi*(m+1), gj*(m+1));   // Do the swap
            }
        }
    }

    newStage();
    //finalize(); // Set all centres correctly
#ifdef SAVE_DATA
    fclose(fsol);
#endif // SAVE_DATA
}
