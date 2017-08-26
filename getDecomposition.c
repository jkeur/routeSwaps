/** Functions to get an optimal cycle decomposition of the move graph G'
 * belonging to the problem of routing qubits from distribution x to y
 * in the fully connected star graph G = (V, E).
 *
 * Created by J. Keur
 * 170818
 */

#define EOC             (unsigned)-1    // End Of Cycle
//#define OK              (unsigned)-1    //
#define LINE_LEN        50              // Length of a line from . file
// Console text colors
#define RED             12              // Red text color
#define NORMALC         7               // Usual text color
#define COLOR_TEXT      SetConsoleTextAttribute(hConsole, RED)      // Red text
#define COLOR_NUM       SetConsoleTextAttribute(hConsole, 8)        // Color the swapped numbers
#define NORMAL_TEXT     SetConsoleTextAttribute(hConsole, NORMALC); // Normal text color
// Define centre usage
#define USE_CENTRE      1       // Centre can be used now to swap with a centre
#define USE_LEAF        2       // Centre can be used now to swap with a leaf
#define IGNORE_C        4       // Ignore centre in applying the rules; use it only in finalizing the state
#define NOT_USABLE      (IGNORE_C | (0 << IGNORE_C))    // Centre is not usable now
#define BEING_USED      (IGNORE_C | (1 << IGNORE_C))    // The number of the centre is being swapped
#define CORRECT         (IGNORE_C | (2 << IGNORE_C))    // Swap centre with leaf if necessary
#define SORTED          (IGNORE_C | (4 << IGNORE_C))    // Centre with leafs OK

#include <limits.h>
#include <math.h>

// Function prototypes
extern char swap(const unsigned i, const unsigned j);

static void allocMem();
static char finalize();
static unsigned add(unsigned *v, unsigned *len, const unsigned num);
void setW();
void printW();

unsigned k, m, n;                       // #centres, #leafs/centre, #nodes
unsigned Ns, Nsb;
unsigned depth;
unsigned *x0, *x, *y;                   // Input vector x0, state vector x, output vector y
unsigned **W, **Wc;                     // Move matrix
unsigned *P;                            // Node cover
unsigned *cycle, *ndist, *np;           // Store cycle, distance to node, #paths
unsigned dbg;                           // dbg = 1 to debug the program, otherwise 0
char     *c2use;                        // centres to use (in current stage)
HANDLE   hConsole;

/** Initialize variables
 * 170410 Created
 */
static void init()
{
    getParams();    // Get k, m
    initParams();   // Set n
    allocMem();

    srand(time(NULL));
    hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
}

/** Set x to initial state x0
 * 170530 Created
 */
void setX()
{
    memcpy(x, x0, n * sizeof(unsigned));    // Set the state equal to the initial state
}

static void allocMem()
{
    unsigned i;

    // Allocate memory
    x0 = (unsigned*)malloc(n * sizeof(unsigned));
    if (x0 == NULL)
    {
        puts("Error allocating mem x0");
        exit(EXIT_FAILURE);
    }
    x = (unsigned*)malloc(n * sizeof(unsigned));
    if (x == NULL)
    {
        puts("Error allocating mem x");
        exit(EXIT_FAILURE);
    }
    y = (unsigned*)malloc(n * sizeof(unsigned));
    if (y == NULL)
    {
        puts("Error allocating mem y");
        exit(EXIT_FAILURE);
    }
    W = (unsigned**)malloc(k * sizeof(unsigned*));
    if (W == NULL)
    {
        puts("Error allocating mem W");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < k; i++)
    {
        W[i] = (unsigned*)malloc(k * sizeof(unsigned));
        if (W[i] == NULL)
        {
            puts("Error allocating mem W[]");
            exit(EXIT_FAILURE);
        }
    }
    Wc = (unsigned**)malloc(k * sizeof(unsigned*));
    if (Wc == NULL)
    {
        puts("Error allocating mem Wc");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < k; i++)
    {
        Wc[i] = (unsigned*)malloc(k * sizeof(unsigned));
        if (Wc[i] == NULL)
        {
            puts("Error allocating mem Wc[]");
            exit(EXIT_FAILURE);
        }
    }
    ndist = (unsigned*)malloc(k * sizeof(unsigned));
    if (ndist == NULL)
        exit(EXIT_FAILURE);
    np    = (unsigned*)malloc(k * sizeof(unsigned));
    if (np == NULL)
        exit(EXIT_FAILURE);
    cycle = (unsigned*)malloc(k * sizeof(unsigned));
    if (cycle == NULL)
        exit(EXIT_FAILURE);
    P     = (unsigned*)malloc(k * sizeof(unsigned));
    if (P == NULL)
    {
        puts("Error allocating P");
        exit(EXIT_FAILURE);
    }
    c2use = (char*)malloc(k * sizeof(char));
    if (c2use == NULL)
        exit(EXIT_FAILURE);
}

/** Reset the variables for this round (depth, Ns, Nsb)
 * 170520 Created
 */
void newRound()
{
    unsigned i;

    // Initialize round specific variables
    depth = 0;          //
    Ns    = 0;          // No swaps done
    Nsb   = 0;          // No expensive swaps done with cost b
    dbg   = 0;

    memset(c2use, 1, k * sizeof(char));
    for (i = 0; i < k; i++)
    {
        //memset(simpl[i], 0, k * sizeof(unsigned));
        c2use[i] = USE_CENTRE | USE_LEAF;
    }

    setX();
    setW();
}

/** Initialize some and update other variables to begin a new stage
 * 1705 Created
 */
static void newStage()
{
    unsigned gj;

    // Reset the state indication; each node can be used in any swap
    for (gj = 0; gj < k; gj++)
    {
        if (c2use[gj] < CORRECT)
            c2use[gj] = USE_CENTRE | USE_LEAF;
        else if (c2use[gj] == (BEING_USED | CORRECT))
            c2use[gj] = CORRECT;
    }
    depth++;
    //swapDone = 0;
    if (dbg)
    {
        printf("\n> d%*u: ", Nd, depth);
    }

    setW();
    finalize();
}

/** For each group, check if all numbers are there and if the centre node
 * does not have the correct number, do a swap.
 * It returns a non-zero value if all numbers are sorted, otherwise it returns 0.
 * 170510 Created
 */
static char finalize()
{
    unsigned i, l;
    char busy = 0;

    for (i = 0; i < k; i++)
    {
        if (c2use[i] == CORRECT)    // If the right number should be brought from a leaf to the centre
        {
            for (l = 1; l <= m; l++)
                if (x[i*(m+1) + l] == i*(m+1))
                {
                    if (swap(i*(m+1), i*(m+1) + l)) // If the swap with this leaf succeeded
                        c2use[i] = SORTED; // This node group is sorted now
                    break;
                }
        }

        // Check if all numbers are in the correct place
        if (c2use[i] != SORTED)     // If not all numbers in group Gi are correct
            busy = 1;
    }
    //setW();     // Recover W

    return busy;
}

/** Check if nodes gi, gj are in >= 1 cycle
 * 170607 Created
 */
static unsigned handleCycle(const unsigned gi, const unsigned gj)
{
    unsigned cnt;

    if (W[gi][gj])
    {
        cnt = (unsigned)fmin((float)W[gi][gj], (float)W[gj][gi]);   // So many efficient swaps are possible over edge (i,j)
        W[gi][gj] -= cnt;
        W[gj][gi] -= cnt;
        return cnt;
    }
    return 0;
}

/** Delete 2-cycles from W & count #2-cycles
 * 170701 Created
 */
unsigned del2cycles()
{
    unsigned gi, gj, cnt = 0;

    for (gi = 0; gi < k; gi++)
        for (gj = gi + 1; gj < k; gj++)
            cnt += handleCycle(gi, gj); // Handle cycles of length 2: edge (i,j)
    return cnt;                         // Return #2-cycles
}

/** Check if W empty on non-diagonal entries
 * 170803 Created
 */
char emptyGraph()
{
    unsigned i, j;

    for (i = 0; i < k; i++)
        for (j = 0; j < k; j++)
        {
            if (W[i][j] > m + 1)
            {
                printf("! W(%u,%u) > %2u\n", i+1, j+1, m+1);
                printW();
                getchar();
            }
            if (i != j)
            {
                if (W[i][j])
                    return 0;
            }
        }
    return 1;
}

/** Get the minimum cycle length
 * 170610 Created
 * 170612 Works correctly
 */
static unsigned getMinCycleLen()
{
    unsigned gs, gi, gj, len, done, Niter = 0;  // Start node s
    unsigned edgesOut, lmin = UINT_MAX;

    if (W[0][0] > m + 1)                        // If W is not set yet
    {
        puts("gMCL:\tW set");
        setW();
    }
    for (gs = 0; gs < k; gs++)
    {
        memset(ndist, 0, k * sizeof(unsigned));
        len = 0;
        done = 0;
        while (!done)
        {
            if (Niter == k*k)
            {
                printW();
                puts("! gMCL:\tCheck this!");
                getchar();
            }
            Niter++;
            for (gi = 0; gi < k; gi++)          // Keep walking until a cycle has been walked
            {
                if (len == 0)
                    gi = gs;
                if (ndist[gi] == len)
                {
                    edgesOut = 0;
                    for (gj = 0; gj < k; gj++)  // Look for nodes j which can be reached from node i
                    {
                        if ( (gj == gi) || (W[gi][gj] == 0) || ndist[gj] )
                            continue;
                        edgesOut = 1;           // Node i has >= 1 outgoing edges
                        // There exists >= 1 path from node i->j && the distance to node j is 0
                        if ( (len >= 1) && (gj == gs) ) // If a cycle has been walked
                        {
                            if (len + 1 < lmin)
                                lmin = len + 1; // Track minimum cycle length
                            gi = k;         // Break outer loop gi; search cycles for the next node s
                            done = 1;       // Min length of cycle found for node s; continue with the next node s
                            break;
                        }
                        ndist[gj] = len + 1;// Set distance from node s->j
                    }
                    if (len == 0)
                    {
                        if (!edgesOut)      // If node i does not have outgoing edges
                            done = 1;
                        break;
                    }
                }
            } // End for gi
            len++;                          // Increase path length
        }
    }

    return lmin;
}

/** Return length of shortest path from pi --> pj iff there are <= Wij shortest paths, otherwise return 0
 * 170802 Created
 */
unsigned wShortestPaths(const unsigned pi, const unsigned pj)
{
    unsigned len = 0, plen;
    unsigned p1, p2;

    memset(ndist, 0, k * sizeof(unsigned));
    memset(np, 0, k * sizeof(unsigned));
    memset(cycle, UINT_MAX, k * sizeof(unsigned));

    while (np[pj] == 0)                     // While no shortest path is found
    {
        for (p1 = 0; p1 < k; p1++)              // Keep walking until a cycle has been walked
        {
            if (len == 0)
            {
                p1 = pi;
                np[pi] = W[pj][pi];
            }
            if (ndist[p1] == len)
            {
                for (p2 = 0; p2 < k; p2++)      // Look for nodes p2 which can be reached from node p1
                {
                    if ( (p2 == p1) || (W[p1][p2] == 0) )
                        continue;
                    // There exists >= 1 path from node p1 -> p2
                    if ( (len >= 1) && (p2 == pj) ) // If a shortest path has been found
                    {
                        np[p2] += (unsigned)fminf((float)np[p1], (float)W[p1][p2]);
                    }
                    else if ( (ndist[p2] == 0) || (ndist[p2] == len + 1) )
                    {
                        np[p2] += (unsigned)fminf((float)np[p1], (float)W[p1][p2]);
                        ndist[p2] = len + 1;    // Set distance from node pi -> p2
                    }
                }
                if (len == 0)
                    break;
            }
        } // End for p1
        if (len >= n)
        {
            printW();
            printf("! oSP: len: %u/%u\n", len, n);
            getchar();
        }
        len++;                                  // Increase path length
    } // End while npj==0

    // Backtrack a shortest cycle & store the cycle
    cycle[0] = pi;
    cycle[len] = pj;
    p2 = pj;
    plen = len;
    len--;
    for (p1 = 0; len && (p1 < k) && (cycle[len] == EOC); )
    {
        if ( (p1 != p2) && W[p1][p2] && (ndist[p1] == len) )
        {
            cycle[len--] = p1;
            p2 = p1;
            p1 = 0;
        }
        else
            p1++;
    }

    if (np[pj] <= np[pi])           // If <= Wij shortest path from pi->pj
        return plen;                // Return length of path
    else
        return 0;                   // Multiple shortest paths from pi->pj
}

/** Get #outgoing edges of centre gi
 * NOTE: W should be set in advance
 * 170607 Created
 */
static unsigned getEDegOut(const unsigned gi)
{
    unsigned gj, deg = 0;

    for (gj = 0; gj < k; gj++)
        if ( (gj != gi) && (W[gi][gj]) )
            deg++;   // #outgoing edges

    return deg;
}

/** Delete a cycle under certain conditions (cond)
 * cond = 0: Remove a shortest cycle if the shortest path from gs to g2 is unique & if gs has 1 outgoing edge
 * cond = 1: Remove a cycle if the shortest path from gs to g2 is unique & if gs has 1 outgoing edge
 * cond = 2: Remove a shortest cycle if the shortest path from gs to g2 is unique
 * cond = 3: Remove a cycle if the shortest path from gs to g2 is unique
 * cond = 4: Remove a shortest cycle if the 'flow' f from g2 to gs satisfies f <= w(gs,g2) & if gs has 1 outgoing edge
 * cond = 5: Remove a shortest cycle if the 'flow' f from g2 to gs satisfies f <= w(gs,g2)
 * cond = 6: Remove a shortest cycle
 * cond = 7: Remove a cycle
 *
 * Return 1 if some cycle is deleted, otherwise return 0 (no cycle deleted)
 * 170612 Created
 * 170613 Worked well
 * 170616 Extended with different conditions determined by cond
 */
static unsigned delCycle(const unsigned gs, const unsigned lmin, const char cond)
{
    unsigned g2, gi, gj, len, cnt, Ci;

    for (g2 = 0; g2 < k; g2++)                  // Keep walking until a cycle has been walked
    {
        if ( (g2 == gs) || (W[gs][g2] == 0) )
            continue;
        memset(ndist, 0, k * sizeof(unsigned));
        memset(np, 0, k * sizeof(unsigned));
        ndist[g2] = 1;                          // Check edge (gs,g2)
        len = 1;
        cnt = 0;                                // Reset cycle counter
        while (np[gs] == 0)                     // While no cycle has been found
        {
            for (gi = 0; gi < k; gi++)          // Keep walking until a cycle has been walked
            {
                if ( (len == 1) && (gi != gs) )
                {
                    gi = g2;
                    np[g2] = 1;                 // Initialize #paths from node s->.. to 1
                }
                if (ndist[gi] == len)
                {
                    for (gj = 0; gj < k; gj++)  // Look for nodes j which can be reached from node i
                    {
                        if ( (gj == gi) || (W[gi][gj] == 0) )
                            continue;
                        // There exists >= 1 path from node i->j
                        if ( (len >= 1) && (gj == gs) ) // If a cycle has been walked
                        {
                            if (cond == 6)
                                np[gj]++;       // Node j can be reached over +1 edge
                            else
                                np[gj] += (unsigned)fminf((float)np[gi], (float)W[gi][gj]);
                            cnt++;              // Count #cycles, not weighted
                        }
                        else if ( (ndist[gj] == 0) || (ndist[gj] == len + 1) )// gj != gs should hold
                        {
                            if (ndist[gj] == len + 1)
                                cnt++;          // Increase #paths
                            ndist[gj] = len + 1;// Set distance from node s->j
                            if (cond == 6)
                                np[gj]++;       // Node j can be reached over +1 edge
                            else
                                np[gj] += (unsigned)fminf((float)np[gi], (float)W[gi][gj]);
                        }
                    }
                    if (len == 1)
                        break;
                }
            } // End for gi
            if (len >= n)
            {
                printW();
                printf("! dC: len: %u/%u\n", len, n);
                getchar();
            }
            len++;                              // Increase path lengths
        } // End while cnt==0

        if ( ( (cond == 0) && (cnt == 1) && (getEDegOut(gs) == 1) && (len == lmin) )
                || ( (cond == 1) && (cnt == 1) && (getEDegOut(gs) == 1) ) // Cond. 0 and 1: remove the shortest cycle containing edge (gs,g2) having a unique path from g2 to gs
                || ( (cond == 2) && (cnt == 1) && (len == lmin) ) // If there is exactly 1 cycle found with length lmin
                || ( (cond == 3) && (cnt == 1) )
                || ( (cond == 4) && (len == lmin) && (np[gs] <= W[gs][g2]) && (getEDegOut(gs) == 1) )
                || ( (cond == 5) && (len == lmin) && (np[gs] <= W[gs][g2]) )
                || ( (cond == 6) && (len == lmin) )
                || (cond == 7) )                // OR if some cycle should absolutely be removed
        {
            if (dbg)
                printf("> Del cycle(%*u,%*u) %u\n", Nd, gs+1, Nd, g2+1, len);

            // Remove 1 cycle by backtracking
            gj = gs;
            cnt = len;                          // Temp. store len
            len--;
            Ci = 0;
            for (gi = 0; gi < k; gi++)
            {
                if ( (gi != gj) && W[gi][gj] && (ndist[gi] == len) ) // If there exists a path from node i->j
                {
                    if (dbg)
                        printf("%*u <- ", Nd, gj+1);
                    cycle[Ci++] = gj;           // Store nodes in cycle (in opposite direction)
                    W[gi][gj]--;                // Remove edge from cycle
                    len--;
                    if (len == 0)
                        break;
                    gj = gi;                    // Backtrack further
                    gi = -1;                    // Search the previous node in the cycle
                }
            }
            if (dbg)
                printf("%*u\n", Nd, gi+1);
            len = cnt;                          // Restore cycle length
            cycle[Ci++] = gi;                   // Store nodes in cycle (in opposite direction)
            W[gs][gi]--;                        // Remove last edge from cycle
            return 1;                           // +1 cycle found & deleted
        }
    }
    return 0;                                   // No cycle deleted
}

/** Delete cycles of length [len] from W
 * 170803 Created
 */
unsigned delCycles(const unsigned len)
{
    unsigned i, gi, pj, gs, ci, cnt;
    unsigned Nc = 0;                        // #cycles removed

    if (len == 2)
    {
        // Delete 2-cycles
        for (gi = 0; gi < k; gi++)
            for (pj = gi + 1; pj < k; pj++)
                if ((ci = handleCycle(gi, pj))) // If there is >= 1 2-cycle => remove from W
                    Nc += ci;               // Count #cycles removed

        return Nc;                          // Ready
    }

    // Delete k-cycles, k >= 3
    memset(np, 0, k * sizeof(unsigned));    // Set np[i] if node i visited
    memset(cycle, UINT_MAX, k * sizeof(unsigned));
    for (gs = 0; gs < k; gs++)
    {
        np[gs] = 1;                         // Node gs is visited now
        cycle[0] = gs;                      // Start walking from here
        cycle[1] = EOC;
        cnt = 0;

        for (i = 1; i; )
        {
            gi = cycle[i] + 1;              // Search further from this node

            for (; gi < k; gi++)            // Search next step
            {
                if (i == len)
                {
                    if (W[cycle[i-1]][gs])
                        gi = gs;
                    else
                        break;
                }
                if ( (cycle[i-1] == gi) || ((i < len) && (cycle[i] == gi)) || !W[cycle[i-1]][gi] )
                    continue;
                // Take a step further in the walk
                np[gi] = 1;                 // Node i visited
                if (Wc[cycle[i-1]][gi])
                {
                    cnt++;                  // Count #edges (i,j) having a unique path from j->i
                }
                if ((ci = add(cycle, &i, gi)) != OK)  // If cycle walked (index i is increased by 1)
                {
                    if (i < k)
                        cycle[i] = EOC;     // Mark end of cycle
                    if ( (i - ci == len) && (cnt >= len - 1) ) // If cycle has length len && if it can be removed
                    {
                        // OK, cycle can be deleted; do it
                        for (pj = 1; pj < len; pj++)
                        {
                            W[cycle[pj-1]][cycle[pj]]--;
                            if (Wc[cycle[pj-1]][cycle[pj]])
                                Wc[cycle[pj-1]][cycle[pj]]--;
                        }
                        W[cycle[pj-1]][cycle[0]]--;
                        if (Wc[cycle[pj-1]][cycle[0]])
                            Wc[cycle[pj-1]][cycle[0]]--;
                        Nc++;                       // +1 cycle removed
                        // Avoid walking edges more often than possible; restart with gs = 0
                        gs = -1;
                        i = 1;              // Break loop i
                        break;
                    }
                    // Take >= 1 step back
                    i--;
                    if (Wc[cycle[i-1]][cycle[i]])
                    {
                        cnt--;
                    }
                    cycle[i] = EOC;
                    break;
                }
                else                            // If no cycle walked (index i is increased by 1)
                {
                    if (i == len + 1)
                    {
                        // Take 1 step back
                        i--;
                        if (Wc[cycle[i-1]][cycle[i]])
                        {
                            cnt--;
                        }
                        cycle[i] = EOC;
                        break;
                    }
                    gi = -1;                    // Start searching a next step
                }
            } // End for gi

            // Take 1 step back
            i--;
            if (i)
            {
                if (Wc[cycle[i-1]][cycle[i]])
                {
                    cnt--;
                }
            }
        } // End for i
    }

    return Nc;                          // #cycles removed
}

/** If node has 1 in-neighbour => simplify W
 * 170807 Created
 */
char ruleB()
{
    unsigned pi, pj, pin, cnt;
    char applied = 0;                       // Rule applied (1) or not (0)

    for (pi = 0; pi < k; pi++)
    {
        cnt = 0;
        for (pj = 0; pj < k; pj++)
        {
            if ( (pj != pi) && W[pj][pi])
            {
                cnt++;
                pin = pj;
            }
        }
        if (cnt == 1)                       // If in-degree = out-degree = 1 => replace adjacent edges
        {
            for (pj = 0; pj < k; pj++)
            {
                if ( (pj != pi) && (pj != pin) && W[pi][pj])   // For each outgoing edge (i,j)
                {
                    W[pin][pi]--;           // Remove edge (in,i)
                    W[pi][pj]--;            // Remove edge (i,j)
                    W[pin][pj]++;           // Add edge (in,j)
                    printf("B(%2u->%2u->%2u)\n", pin+1, pi+1, pj+1);
                    applied = 1;
                }
            }
        }
    }

    return applied;
}

/** If node has 1 out-neighbour => simplify W
 * 170807 Created
 */
char ruleC()
{
    unsigned pi, pj, pout, cnt;
    char applied = 0;                       // Rule applied (1) or not (0)

    for (pi = 0; pi < k; pi++)
    {
        cnt = 0;
        for (pj = 0; pj < k; pj++)
        {
            if ( (pj != pi) && W[pi][pj])
            {
                cnt++;
                pout = pj;
            }
        }
        if (cnt == 1)                       // If in-degree = out-degree = 1 => replace adjacent edges
        {
            for (pj = 0; pj < k; pj++)
            {
                if ( (pj != pi) && (pj != pout) && W[pj][pi])   // For each incoming edge (j,i)
                {
                    W[pj][pi]--;            // Remove edge (j,i)
                    W[pi][pout]--;          // Remove edge (i,out)
                    W[pj][pout]++;          // Add edge (j,out)
                    printf("C(%2u->%2u->%2u)\n", pj+1, pi+1, pout+1);
                    applied = 1;
                }
            }
        }
    }

    return applied;
}


/** Get an optimal cycle decomposition
 * 170802 Created
 */
unsigned getDecomp()
{
    unsigned pi, pj, beta = n, lmin, cnt;
    unsigned Niter;

    for (pi = 0; pi < k; pi++)
        beta -= W[pi][pi];              // These qubits don't have to be moved
    beta -= del2cycles();               // Delete 2-cycles and count them
    for (pi = 0; pi < k; pi++)
        memset(Wc[pi], 0, k * sizeof(unsigned));

    while (!emptyGraph())
    {
        lmin = getMinCycleLen();

        for (Niter = 0; (Niter <= 2) && (getMinCycleLen() == lmin); Niter++)
        {
            for (pi = 0; pi < k; pi++)
            {
                for (pj = 0; pj < k; pj++)
                {
                    if ( (pj != pi) && W[pi][pj] && wShortestPaths(pj, pi) )
                    {
                        Wc[pi][pj] = 1;
                    }
                }
            }
            cnt = delCycles(lmin);
            if (cnt)
            {
                beta -= cnt;
                Niter = -1;
            }
            else if (getMinCycleLen() != lmin)
                break;
            if (Niter == 1)
            {
                if (ruleB())
                    break;
                if (ruleC())
                    break;
            }
            else if (Niter == 2)
            {
                // 'Randomly' remove one of the cycles
                for (pi = 0; (pi < k) && !delCycle(pi, lmin, 2); pi++);
                if (pi < k)         // If success
                    beta--;
                else                // If no success
                {
                    for (pi = 0; (pi < k) && !delCycle(pi, lmin, 5); pi++);
                    if (pi < k)     // If success
                        beta--;
                    else
                    {
                        for (pi = 0; (pi < k) && !delCycle(pi, lmin, 6); pi++);
                        if (pi < k)     // If success
                            beta--;
                        else
                        {
                            printW();
                            puts("No success!");
                            getchar();
                        }
                    }
                }
                break;
            }
        } // End for Niter
    }
    setW();     // Recover W

    return beta;
}

/** Add element to vector
 * Returns OK if added, element index if num is in v already
 * 170509 Created
 */
static unsigned add(unsigned *v, unsigned *len, const unsigned num)
{
    unsigned i;

    for (i = 0; i < *len; i++)
    {
        if (v[i] == num)
            return i;               // Element is in the vector v; do not add it
    }
    v[(*len)++] = num;
    return OK;
}

/** Get star destination of xj
 * 170711 Created
 */
unsigned getDestStar(const unsigned j)
{
    const int posYj = inv(y, x[j]);
    if (posYj == -1)
    {
        puts("! posYj = -1");
        exit(EXIT_FAILURE);
    }
    return posYj/(m+1);
}

/** Test if number xj at node j has a centre destination or not.
 * 170711 Created
 */
char destIsCentre(const unsigned j)
{
    return inv(y, x[j]) % (m+1) == 0;
}

/** Set move matrix W and centre move matrix Wc
 * 170411 Created
 */
void setW()
{
    unsigned i;
    int di;

    if (x[0] > n)
    {
        memcpy(x, x0, n * sizeof(unsigned));
        puts("! setW: x was not set");
        getchar();
    }
    for (i = 0; i < k; i++)
    {
        memset(W[i], 0, k * sizeof(unsigned));
        memset(Wc[i], 0, k * sizeof(unsigned));
    }
    for (i = 0; i < n; i++)
    {
        if (x[i] == 0)
            continue;
        di = getDestStar(i);                // Get destination
        if (di == -1)
        {
            puts("! di = -1");
            exit(EXIT_FAILURE);
        }
        else if (di > n)
        {
            printf("! di: %u\n", di);
            exit(EXIT_FAILURE);
        }
        W[i/(m+1)][di]++;
        if (destIsCentre(i))                // If xi has a centre destination
            Wc[i/(m+1)][di]++;
    }
}

/** Print move matrix W
 * 170412 Created
 */
void printW()
{
    unsigned i, j;

    printf("W =\n%*c | ", Nd, 'c');
    for (i = 1; i <= k; i++)
        printf("%*u ", Nd, i);
    puts("");
    for (i = 0; i < k; i++)
    {
        printf("%*u | ", Nd, i+1);
        for (j = 0; j < k; j++)
            printf("%*u ", Nd, W[i][j]);
        puts("");
    }

    printf("Wc =\n%*c | ", Nd, 'c');
    for (i = 1; i <= k; i++)
        printf("%*u ", Nd, i);
    puts("");
    for (i = 0; i < k; i++)
    {
        printf("%*u | ", Nd, i+1);
        for (j = 0; j < k; j++)
            printf("%*u ", Nd, Wc[i][j]);
        puts("");
    }
}

/** Calculate #moves: a = 1^T*W*1 - sum_i W_ii = n - sum(.)
 * 170523 Created
 */
unsigned getNmoves()
{
    unsigned gi, alpha = n;

    for (gi = 0; gi < k; gi++)
        alpha -= W[gi][gi];

    return alpha;
}
