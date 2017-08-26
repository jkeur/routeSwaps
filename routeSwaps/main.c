/** A program to route swaps in Fully Connected Star graphs
 *
 * Created by J. Keur
 * 170816
 *
 * Problem: Sort input sequence x to permuted output sequence y using a fully connected star graph containing
 * k centres and m leaves per centre, where node j = i(m+1)+1 (with i in [0,k-1], j in [1,k(m+1)])
 * is a centre node and the others are leaves.
 */

#define REPEAT  500             // Repeat the protocol ... times
//#define LOAD_P

#include <stdio.h>
#include <stdlib.h>
#include <windows.h>            // For SYSTEMTIME
#include <time.h>
#include <sys/time.h>
#include "../handleVariables.c"
#include "../getDecomposition.c"
#include "../routeSimple.c"

extern unsigned *x0, *y;        // Input vector, output vector
extern unsigned depth, Ns, Nsb;
extern unsigned dbg;            // dbg = 1 to debug the program, otherwise 0

unsigned        Nd;
unsigned        start, finish;
SYSTEMTIME      st;
FILETIME        ft;
ULARGE_INTEGER  ui;

// Function prototypes
void tic();
unsigned toc();
static void loadP(const char *fname);
static void solveByMyAlg(const char *fname);
void load(const char *fname);
void save(char *fname);

int main()
{
    static unsigned DtT, NstT, NsbtT;   // For trivial alg.: Total #stages, total #swaps, total #expensive swaps
    static unsigned DtM, NstM, NsbtM;   // For my alg.:      Total #stages, total #swaps, total #expensive swaps
    static unsigned msecT, msecM;
    unsigned i, a, b, b2, cnt0;
    unsigned nonOpt;                    // #non-optimal solutions
    char in;
    char fname[LINE_LEN] = "p";

    init();
    dbg   = 0;
    DtT   = 0;
    NstT  = 0;
    NsbtT = 0;
    DtM   = 0;
    NstM  = 0;
    NsbtM = 0;
    a     = 0;
    b     = 0;
    msecT = 0;
    msecM = 0;
    cnt0  = 0;  // Count #times alg. 1 is better than 2
    nonOpt = 0;
    for (i = 0; i < n; i++)
        y[i] = i + 1;

#ifdef LOAD_P
    sprintf(fname, "p");
    loadP(fname);
    dbg = 1;
    b2     = getDecomp();
    b      = pealW();
    dbg = 0;
    printf("b: %u | b2: %u\n", b, b2);

    //setX();
    printW();
    dbgW();
    del2cycles();
    for (i = 3; i <= 4; i++)
    {
        listCycles(i);
        printC();
    }
    printNcLen();
    return 0;
#endif // LOAD_P

    for (i = 0; i < REPEAT; i++)
    {
        setRandom(x0);
#ifdef PRINT_STATE
        puts("> BEGIN state:");
        printState();
        setW();
        setD();
        printW();
        printD();
#endif // PRINT_STATE

        // Run trivial algorithm & set statistics
        tic();
        routeSimple();                   // RUN trivial algorithm routeSimple
        msecT += toc();
        DtT   += depth;
        NstT  += Ns;
        NsbtT += Nsb;
#ifdef PRINT_STATE
        printf("> FINAL state trivial alg.:");
        printState();
#endif // PRINT_STATE
        // Run my algorithm & set statistics
        setX();
        setW();
        tic();
        b2 = getDecomp();
        msecM += toc();
        b     += b2;
        a     += getNmoves();
        DtM   += depth;
        NstM  += Ns;
        NsbtM += Nsb;
#ifdef PRINT_STATE
        printf("> FINAL state my alg.:");
        printState();
        printW();
#endif // PRINT_STATE
    }
    printf("\n%u/%u\n", cnt0, REPEAT);
    NstM  = a;
    NsbtM = b;

    // Print statistics
    printf("\nAVG %u\tSimple\tMy\tBetter\tWon\tOpt\tMy-Opt\n", REPEAT);
    printf("d    \t%.1f\t%.1f\t%c\t%.1f\n", (float)DtT/REPEAT, (float)DtM/REPEAT, DtM < DtT ? 'Y' : ' ', ((float)DtT - (float)DtM)/REPEAT);
    printf("#s(a)\t%.1f\t%.1f\t%c\t%.1f\t%.1f\n", (float)(NstT-NsbtT)/REPEAT, (float)(NstM-NsbtM)/REPEAT, (NstM-NsbtM) < (NstT-NsbtT) ? 'Y' : ' ', ((float)(NstT-NsbtT)-(float)(NstM-NsbtM))/REPEAT, (float)a/REPEAT);
    COLOR_TEXT;
    printf("#s(b)\t%.1f\t%.1f\t%c\t%.1f\t%.1f\t%.1f\n", (float)NsbtT/REPEAT, (float)NsbtM/REPEAT, NsbtM < NsbtT ? 'Y' : ' ', (float)(NsbtT - NsbtM)/REPEAT, (float)b/REPEAT, ((float)NsbtM-(float)b)/REPEAT);
    NORMAL_TEXT;
    printf("#s   \t%.1f\t%.1f\t%c\t%.1f\n", (float)NstT/REPEAT, (float)NstM/REPEAT, NstM < NstT ? 'Y' : ' ', ((float)NstT - (float)NstM)/REPEAT);
    printf("M    \t \t \t \t \t%.1f\n", (float)NstM/REPEAT);
    printf("time \t%.2f\t%.2f\n", (float)msecT/REPEAT, (float)msecM/REPEAT);
    printf("Non-opt\t%2u\n", nonOpt);

    while ((in = getchar()) != 'c')
    {
        switch (in)
        {
        case 'e':           // Solve the problem using the trivial algorithm
            tic();
            routeSimple();
            printf("\trouteSimple\tOpt\nd\t%u\n#s(a)\t%u\n#s(b)\t%u\n#s\t%u\ntime\t%u ms\n", depth, Ns-Nsb, Nsb, Ns, toc());
            break;

        case 'h':
            printf("Press one of the following keys\nc: Close the solver\ne: Use the trivial solving algorithm\nh: Show help information\nm: Use my advanced algorithm\ns: Save the problem\n");
            break;

        case 'm':
            solveByMyAlg(fname);
            break;

        case 's':
            save(fname);    // Save fully connected star graph problem
            break;

        case 'l':
            sprintf(fname, "p4");
            loadP(fname);   // Load fully connected star graph problem
        }
    }

    return 0;
}

/** Reset timer. The timer value is obtained by calling toc().
 * 170523 Created
 */
void tic()
{
    GetSystemTime(&st);
    SystemTimeToFileTime(&st, &ft); // converts to file time format
    ui.LowPart  = ft.dwLowDateTime;
    ui.HighPart = ft.dwHighDateTime;
    start = ui.QuadPart;
}

/** Register time elapse after the last call of tic()
 * 170523 Created
 */
unsigned toc()
{
    GetSystemTime(&st);
    SystemTimeToFileTime(&st, &ft); // converts to file time format
    ui.LowPart  = ft.dwLowDateTime;
    ui.HighPart = ft.dwHighDateTime;
    finish = ui.QuadPart;
    return (unsigned)(((float)(finish-start))/10000);
}

/** Load problem
 * 170624 Created
 */
static void loadP(const char *fname)
{
    load(fname);        // Load fully connected star graph problem
    setX();
}

/** Solve the problem with my solving algorithm
 * 170624 Created
 */
static void solveByMyAlg(const char *fname)
{
    unsigned t;
    tic();
    t = toc();
    printf("\tRouteSwaps\nd\t%u\n#s(a)\t%u\n", depth, Ns-Nsb);
    COLOR_TEXT;
    printf("#s(b)\t%u\n", Nsb);
    NORMAL_TEXT;
    printf("#s\t%u\ntime\t%u ms\n", Ns, t);
}

/** Load problem from path
 * 170523 Created
 */
void load(const char *fname)
{
    unsigned i;
    char line[LINE_LEN];
    char path[LINE_LEN];

    sprintf(path, "%s.fcs", fname);     // Set path of the problem
    FILE *fp = fopen(path, "rt");
    if (fp == NULL)
    {
        printf("! Failed to open \"%s\"\n", path);
        return;
    }

    while ( fgets(line, LINE_LEN, fp) && (line[0] != 'p') );    // Find problem description
    sscanf(&line[2], "%u %u", &k, &m);  // Read parameters k, m
    n = k * (m+1);                      // Compute #nodes
    while ( fgets(line, LINE_LEN, fp) && (line[0] != 'x') );    // Find initial assignment x0
    for (i = 0; fgets(line, LINE_LEN, fp); i++)
    {
        sscanf(line, "%u", &x0[i]);     // Read state numbers in [n]
        x0[i]--;                        // Convert numbers to [n]-1
    }

    fclose(fp);
    printf("> Problem loaded from \"%s\"\n", path);
}

/** Save problem as '[fname].fcs'
 * 170523 Created
 */
void save(char *fname)
{
    unsigned i;
    char path[LINE_LEN];
    sprintf(path, "%s.fcs", fname);
    FILE *fp = fopen(path, "w");
    if (fp == NULL)
    {
        printf("! Unable to use the path \"%s\"\n", path);
        return;
    }
    fprintf(fp, "c Generated by sortFuConStar\nc Used %u centres with %u leafs each: %u nodes\np %u %u\nx\n", k, m, n, k, m);
    for (i = 0; i < n; i++)
        fprintf(fp, "%*u 0\n", Nd, x0[i] + 1);  // Save initial input vector x0

    fclose(fp);
    printf("> Problem saved as \"%s\"\n", path);
}
