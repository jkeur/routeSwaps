/** Get input values for sorting networks/sequences
 *
 * Created by J. Keur
 *
 * 170410
 */

#define SET_M   5               // Set m to a predefined value: 5 leaves/star

#include <limits.h>
#include <math.h>

#define OK                  UINT_MAX

// External variables
extern unsigned k, m, n, Nd;    // #stars, #nodes/star, total #nodes, #digits to represent node labels
extern unsigned Nc;             // #clauses, #digits of decimal variable number

// Function prototypes
int inv(const unsigned *y, const unsigned yi);

/** Shuffle vector x
 * 170414 Created
 */
void setRandom(unsigned *x)
{
    unsigned i, v1, v2, valT;

    for (i = 0; i < n; i++)
        x[i] = i + 1;
    for (i = n*n/2; i; i--)
    {
        // Randomly shuffle current node values
        v1 = rand() % n;
        v2 = rand() % n;            // Randomly select node numbers to shuffle value thereof
        valT = x[v1];               // Temporarily store value
        x[v1] = x[v2];
        x[v2] = valT;
    }
}

/** Get input sequence of numbers
 * 170321 Created
 * 170407 Removed fprintf to SAT file
 */
void getVals(unsigned *vals, const unsigned minOne)
{
    unsigned i, ni, val, x;
    char in[10];

    printf("? Give the value of each line (in [1, %u])\n", n);
    x = 1;
    for (ni = 0; ni < n; )
    {
        printf("%2u: ", ni + 1);
        for (i = 0; ((in[i] = getchar()) != ' ') && (in[i] != '\n'); i++);
        in[i] = 0;
        if ((in[0] >= '0') && (in[0] <= '9'))
        {
            val = atoi(in);
            if (minOne)
                val--;
            if (val >= n)
            {
                printf("! %u >= n = %u, I'll make this 0\n", val, n);
                vals[ni] = 0;
            }
            else
                vals[ni] = val;
        }
        else
        {
            // No number entered
            if (in[0] == 0)
                continue;
            switch (in[0])
            {
            case 'i':
            {
                if (x > 1)
                {
                    puts("! Only possible if no input is given before");
                    break;
                }
                // Set input in increasing order: (n-1, n-2, ..., 1, 0)
                for (i = 0; i < n; i++)
                    vals[i] = i;
                printf("> Input values set in increasing order: (1, ..., %u)\n", n);
            }
            return;
            case 'd':
            {
                if (x > 1)
                {
                    puts("! Only possible if no input is given before");
                    break;
                }
                // Set input in decreasing order: (n-1, n-2, ..., 1, 0)
                for (i = 0, val = n - 1; i < n; i++, val = n - i - 1)
                {
                    vals[i] = val;
                }
                printf("> Input values set in decreasing order: (%u, ..., 1)\n", n);
            }
            return;
            case 'r':
            {
                setRandom(vals);
                puts("> Input vales set randomly");
            }
            return;
            case 'h':
            {
                // Ask for help
                puts("> Use the following inputs:\n- A number\n- 'd': Set input values in decreasing order");
            }
            break;
            default:
                printf("! Unknown command: %c (%u)\n", in[0], in[0]);
            }
        }
        ni++;
    }
}

/** Check if values are unique and in [1, n]
 * 170411 Created
 */
unsigned checkVals(const unsigned *vals)
{
    unsigned i;

    for (i = 0; i < n; i++)
        if (inv(vals, i) == -1)
            return i;   // No inverse mapping found for number i
    return OK;
}

/** Get inverse mapping of c(i) = y_i
 * 170410 Created
 */
int inv(const unsigned *y, const unsigned yi)
{
    unsigned i;

    for (i = 0; i < n; i++)
        if (y[i] == yi)
            return i;

    return -1;          // No inverse mapping found
}

/** Print sorted input values
 * 170321 Created
 */
void printVals(unsigned *vals)
{
    unsigned i;

    printf("  n val\n-------\n");
    for (i = 0; i < n; i++)
        printf("%*u %*u\n", Nd, i + 1, Nd, vals[i]);
}


/** Functions to get the parameters k, m
 * and set the parameters n, Nd.
 */

/** Get the parameters k and m
 * 170524 Created
 */
void getParams()
{
    k = 0;
    // Get #centres k
    while ((k == 0) || (k > 100))
    {
        printf("? Give #centres k: ");
        scanf("%u", &k);
    }
#ifdef SET_M
    m = SET_M;
    printf("> Using %u leafs/centre\n", SET_M);
#else
    m = 0;
    // Get #leafs/centre m
    while ((m == 0) || (m > 100))
    {
        printf("? Give #leafs/c m: ");
        scanf("%u", &m);
    }
#endif // SET_M
}

/** Set the parameters n, e, Nd, using k and m
 * 170524 Created
 */
void initParams()
{
    n = k*(m+1);                // #nodes
    Nd = ceil(log10(n+1));      // #digits of node numbers
}
