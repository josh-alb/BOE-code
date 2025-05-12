/* This is a probabilistic meter-finding program. It takes "notelist"
input (a series of statements of the form "Note [ontime] [offtime]
[pitch]", with ontime and offtime in milliseconds and pitch in integer
notation, middle C = 60). It assumes monophonic input.  It derives a
three-level metrical grid. Parameters can be set at the top of the
code below; the values shown are taken from the Essen folksong
database.

Compile it like this: "cc -lm meter16.c -o meter16". Run it like this:
"./meter16 [input-file]"

The output can be controlled with the -v flag. With -v = -2, the
output is just the probability of the rhythmic pattern. With -v = -1,
it's a list of"Note" and "Beat" statements. (A Beat statement has the form
"Beat [time] [level]", indicating the metrical level of each beat. This 
format can be used as input to the Melisma programs, and also to generate 
a note-address list.) With -v = 1, it's a graphic display of the metrical 
structure with the notes.

The workings of the program are explained in _Music and Probability_,
chapter 3.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/***** Parameters *****/

#define UW 0.95       /* P of a note on a level 3 beat ("upper weight") */
#define TW 0.74       /* P of a note on a level 2 beat ("tactus weight") */
#define LW 0.38       /* P of a note on a level 1 beat ("lower weight") */
#define NB 0.01      /* P of a note on a pip that's not a beat at all */

#define LDP 0.78      /* P of the lower level being duple (1-LDP is the P that it's triple) */
#define UDP 0.76     /* P of the upper level being duple (1-UDP is the P that it's triple) */

#define DZP 0.65     /* P that the duple upper level starts at phase zero (i.e. the first L2 beat is an L3 beat). 
			(1-DZP is the P that it	starts at phase 1) */
#define TZP 0.33     /* P that the triple upper level starts at phase 0 */
#define TTP 0.996    /* P that, if the triple upper level doesn't start at phase 0, it
			 starts at phase 2 (i.e. the second L2 beat is a L3 beat); so P of phase 2 is (1-TZP)*TTP. 
			 (1-TZP)(1-TTP) is the P that it starts at phase 1 */

/* The probability of a lower-level interval deviating by x-3 pips from the ideal (halfway or onethirdsway or twothirdsway) point.
(These ideal points are rounded down, i.e. if the ideal point is 3.5 it's rounded to 3.) */
double reg_dist[7] = {.02, .08, .24, .32, .24, .08, .02};
//double reg_dist[7] = {.05, .10, .20, .30, .20, .10, .05};

/* The probability of the first tactus interval being x+9 pips long */
double first_dist[14] = {.10, .20, .30, .23, .13, .03, .006, .002, .001, .0006, .0002, .0001, .00005, .00005};

/***********************/

int v = -1;      // verbosity: -2 = p(surface); -1 = note-beat list; 0 = basic stats; 1 = graphic output;
            // 2 = stats, graphic output, more info

FILE * in_file;
char line[50];
char noteword[10];

struct {
    int ontime;
    int offtime;
    int pitch;
} note[5000];

typedef struct tactus_struct {

    int best_duple;
    double best_duple_score;
    int best_triple1;
    int best_triple2;
    double best_triple_score;
    double score[3];
    int bestk[3];
    double duple_mass;
    double triple_mass;
    double mass[3];

} Tactus;

Tactus tactus[10000][23];

int note_pip[10000];

int final[10000];

int last_pip;

int note_offset;

double duple_score[23][23];
double triple_score[23][23][23];
double grand_total_mass, norm_total_mass;

int numnotes;

double interval_score[14][14];
/* You could just define this as an explicit table of hand-coded values. (We experimented
   with this; see "comments".) Right now, though, it's defined by a function. */

int best_final_j[2][2];
int best_phase[2][2];
double timesig_best_score[2][2];

double make_tables(int, int);
void calculate_duple_scores();
void calculate_triple_scores();
void generate_lower_levels();
void trace_back();
void print_graphic_output();
void print_notebeat_list();
int print_note(int, int);

double sec(int x) {
    //return 1.0;
    return (double)x / 20.0;
}

void calculate_interval_scores() {

    /* Define the probability of tactus interval j given previous interval k */

    int j, k;
    double mass, steepness;

    //return;

    for(k=9; k<23; k++) {
	mass = 0.0;
	for(j=9; j<23; j++) {

	    /* if(abs (k-j) > 3) interval_score[j-9][k-9] = 0.0;
	       else interval_score[j-9][k-9] = reg_dist[ (k-j) + 3]; */

	    //steepness = 10.0 / (double)(k);
	    steepness = 1.0;

	    //interval_score[j-9][k-9] = 1.0 / ((steepness * (double)abs(j-k)) + 1.0);
	    //interval_score[j-9][k-9] += 1.0 / ((double)abs(k-14) + 1.0);
	    //else interval_score[j-9][k-9] = reg_dist[ (k-j) + 3];
	    interval_score[j-9][k-9] = exp(-pow((0.5 * abs(j-k)), 2.0));

	    //if(abs(j-k) > 1) interval_score[j-9][k-9] = 0.0;
	    //Use the above line to prohibit any change of interval of more than 1 or 0 pips 

	    //if(j>14) interval_score[j-9][k-9] = 0.0;
	    //Use the above line to force tactus intervals into a certain range

	    mass += interval_score[j-9][k-9];
	}
	//printf("mass = %6.3f\n", mass);
	for(j=9; j<23; j++) {
	    interval_score[j-9][k-9] = interval_score[j-9][k-9] / mass;
	}
    }

    // The following routine prints out the table of interval scores */
    /*
    for(j=9; j<23; j++) {
	printf("j=%2d: ", j);
	for(k=9; k<23; k++) {
	    printf("%6.2f ", interval_score[j-9][k-9]);
	}
	printf("\n");
	}  */
}
    
int abs(int x) {
    if(x>=0) return x;
    else return 0-x;
}

int main(argc, argv) 
     int argc;
     char *argv[];
{

    int a, n, np, end_time, first_time, i, j;
    int ts, u, best_timesig, best_upper, ph;
    double best_score;

    //verbosity
    v = -1;
    end_time = -1;

    for(np=0; np<10000; np++) note_pip[np]=0;

    a=1;
    in_file=NULL;
    while(1) {
	if(a==argc) break;
	else if(a==argc-1) {
	    if(strcmp(argv[a], "stdin")==0) in_file = stdin;
	    else in_file = fopen(argv[a], "r");
	}
	else if(strcmp(argv[a], "-v")==0) (void) sscanf (argv[a+1], "%d", &v);
	a++;
    }

    if(in_file == NULL) { printf("File not found\n"); exit(1); }

    n = 0;

    while (fgets(line, sizeof(line), in_file) !=NULL) {            
	(void) sscanf (line, "%s", noteword);
	if(line[0] == '\0' || line[0] == '\n') continue;
	if (strcmp (noteword, "Note") == 0) { 
	    (void) sscanf (line, "%s %d %d %d", noteword, &note[n].ontime, &note[n].offtime, &note[n].pitch);
	    n++;
	}
	else if (strcmp (noteword, "End") == 0) { 
	    (void) sscanf (line, "%s %d", noteword, &end_time);
	}
    }

    numnotes = n;

    note_offset = 1050 - note[0].ontime;
    for(n=0; n<numnotes; n++) {
	note[n].ontime += note_offset;
	note[n].offtime += note_offset;
    }
    if(end_time != -1) end_time += note_offset;

    //for(n=0; n<numnotes; n++) printf("Note %d %d %d\n", note[n].ontime, note[n].offtime, note[n].pitch);

    for(n=0; n<numnotes; n++) {
	np = note[n].ontime / 50;
	if(np<10000) note_pip[np]=1;
	else { printf("Error: Note out of time range\n"); exit(1); }
    }
    if(end_time != -1) last_pip = end_time / 50;
    else last_pip = note[numnotes-1].offtime / 50;

    for(i=0; i<=last_pip; i++) {
	for(j=0; j<=last_pip; j++) {
	    tactus[i][j].duple_mass = tactus[i][j].triple_mass = 0.0;
	    for(ph=0; ph<3; ph++) tactus[i][j].mass[ph] = 0.0;
	}
    }

    grand_total_mass=0.0;

    calculate_duple_scores();

    calculate_triple_scores();

    calculate_interval_scores();

    if(v>1) printf("Generating lower levels...\n");
    generate_lower_levels();

    if(v>1) printf("Trying 2/4...\n");
    timesig_best_score[0][0] = make_tables(2, 2);
    if(v>1) printf("Trying 6/8...\n");
    timesig_best_score[1][0] = make_tables(3, 2);
    if(v>1) printf("Trying 3/4...\n");
    timesig_best_score[0][1] = make_tables(2, 3);
    if(v>1) printf("Trying 9/8...\n");
    timesig_best_score[1][1] = make_tables(3, 3);
    /* Rule out 9/8 */
    //timesig_best_score[1][1] = -1000.0;

    if(v>1) printf("2/4 score = %6.3f; 6/8 score = %6.3f; 3/4 score = %6.3f; 9/8 score = %6.3f\n", timesig_best_score[0][0], timesig_best_score[1][0], timesig_best_score[0][1], timesig_best_score[1][1]);

    norm_total_mass = pow(grand_total_mass, (1.0/sec(i)));
    if(v>1) printf("Grand total mass = %6.40f; log mass = %6.6f; per second mass = %6.10f\n\n", grand_total_mass, log(grand_total_mass), norm_total_mass);
    if(v==0) printf("p(surface) = %6.3f\n", log(grand_total_mass));
    if(v==-2) printf("%6.6f\n", log(grand_total_mass));

    best_score = -10000.0;
    for(ts=0; ts<2; ts++) {
	for(u=0; u<2; u++) {
	    if(timesig_best_score[ts][u] > best_score) {
		best_timesig = ts+2;
		best_upper = u+2;
		best_score = timesig_best_score[ts][u];
	    }
	}
    }

    timesig_best_score[best_timesig-2][best_upper-2] = make_tables(best_timesig, best_upper);

    if(v>-1) printf("Best analysis: lower = %d, upper = %d, final (penultimate) phase = %d; score = %6.3f\n", best_timesig, best_upper, best_phase[best_timesig-2][best_upper-2], timesig_best_score[best_timesig-2][best_upper-2]);

    //Here's where you can force it to output a particular meter */
    //best_timesig = 2;
    //best_upper = 3; 
    //timesig_best_score[best_timesig-2][best_upper-2] = make_tables(best_timesig, best_upper);
    //best_phase[0][1] = 0; 

    //printf("max p(structure & surface) = %6.3f\n", timesig_best_score[best_timesig-2][best_upper-2]);
    trace_back(best_timesig, best_upper);
    if(v > 0) print_graphic_output(); 
    else if(v==-1) print_notebeat_list();
}

void calculate_duple_scores() {

    int j, k, halfway_point;

    for(j=9; j<23; j++) {
	for(k=1; k<j; k++) {

	    /*
	    if(j > 18) {
		duple_score[j][k] = -100.0;
		continue;
		} */
	    halfway_point = j / 2;
	    //printf("j = %d, halfway_point = %d\n");
	    if(abs(k-halfway_point) > 3) duple_score[j][k] = -1000.0;
	    else duple_score[j][k] = log(reg_dist[ (k-halfway_point) + 3 ] );
	}
    }
}

	    
void calculate_triple_scores() {

    int j, k1, k2, j2;
    int onethird_point, twothirds_point;
    double mass;

    for(j=9; j<23; j++) {
	mass = 0.0;
	for(k1 = 1; k1 < j-1; k1++) {
	    for(k2 = k1+1; k2 < j; k2++) {
		triple_score[j][k1][k2] = 0.0;
		onethird_point = j / 3;
		if(abs(k1 - onethird_point) > 3) triple_score[j][k1][k2] += -1000.0;
		else triple_score[j][k1][k2] += log(reg_dist[ (onethird_point - k1) + 3 ] );
		twothirds_point = (j+k1) / 2;
		if(abs(k2 - twothirds_point) > 3) triple_score[j][k1][k2] += -1000.0;	    
		else triple_score[j][k1][k2] += log(reg_dist[ (twothirds_point - k2) + 3 ] );
		mass += exp(triple_score[j][k1][k2]);
	    }
	}
	//printf("mass = %6.3f\n", mass);
	for(k1 = 1; k1 < j-1; k1++) {
	    for(k2 = k1+1; k2 < j; k2++) {
		if(mass > 0.001) triple_score[j][k1][k2] = log ((exp (triple_score[j][k1][k2])) / mass);
		else triple_score[j][k1][k2] = -1000.0;
		/* for(j2=0; j2<=j; j2++) {
		    if(j2==k1 || j2==k2) printf("1 ");
		    else printf("0 ");
		    } */
		//printf(" (%d,%d,%d: %6.3f)\n", j, k1, k2, triple_score[j][k1][k2]);
	    }
	}
    }    
}

void generate_lower_levels() {

    int i, j, k, halfway_point, best, k1, k2, bestk1, bestk2, p;
    double best_score;
    double lower_note_score;
    double duple_mass, triple_mass;

    for(i=0; i<=last_pip; i++) {

	for(j=9; j<23; j++) {

	    /* Calculate best lower-level duple beat */

	    best_score = -1000.0;
	    best = -1;

	    for(k=1; k<j; k++) {

		lower_note_score = 0.0;
		for(p=(i-j)+1; p<i; p++) {
		    if(p==(i-j)+k) {
			if(note_pip[p]==1) lower_note_score += log(LW);
			else lower_note_score += log(1.0-LW);
		    }
		    else {
			if(note_pip[p]==1) lower_note_score += log(NB);
			else lower_note_score += log(1.0-NB);
		    }
		}

		//if(i==100) printf("lower_note_score(i=%d,j=%d,k=%d)=%6.3f\n", i, j, k, lower_note_score);
		//if(i==100) printf("(%d, %d): duple_score = %6.3f\n", j, k, duple_score[j][k]);
		if(duple_score[j][k] + lower_note_score > best_score) {
		    best = k;
		    best_score = duple_score[j][k] + lower_note_score;
		}
		tactus[i][j].duple_mass += exp(duple_score[j][k] + lower_note_score);
	    }
	    if(best==-1) {
		printf("Duple Error!\n"); 
		exit(1);
	    }
	    tactus[i][j].best_duple = (i-j) + best;
	    tactus[i][j].best_duple_score = best_score;

	    //if(i==14) printf("Tactus int (%d, %d): best duple k = %d, beat = %d, score = %6.3f\n", i, i-j, best, tactus[i][j].best_duple, tactus[i][j].best_duple_score);

	    /* Calculate best lower-level triple beats */

	    best_score = -10000.0;
	    bestk1 = bestk2 = -1;

	    for(k1=1; k1<j-1; k1++) {
		for(k2=k1+1; k2<j; k2++) {
		    
		    lower_note_score = 0.0;

		    for(p=(i-j+1); p<i; p++) {
			if(p==(i-j+k1) || p==(i-j+k2)) {
			    if(note_pip[p]==1) lower_note_score += log(LW);
			    else lower_note_score += log(1.0-LW);
			}
			else {
			    if(note_pip[p]==1) lower_note_score += log(NB);
			    else lower_note_score += log(1.0-NB);
			}
		    }

		    if( (triple_score[j][k1][k2] + lower_note_score) > best_score) {
			bestk1 = k1;
			bestk2 = k2;
			best_score = triple_score[j][k1][k2] + lower_note_score;
		    }

		    tactus[i][j].triple_mass += exp(triple_score[j][k1][k2] + lower_note_score);
		}
	    }

	    if(bestk1==-1 || bestk2==-1) {
		printf("Triple Error!\n"); 
		exit(1);
	    }

	    tactus[i][j].best_triple1 = (i-j) + bestk1;
	    tactus[i][j].best_triple2 = (i-j) + bestk2;
	    tactus[i][j].best_triple_score = best_score;

	    //printf("Tactus int (%d, %d): best t1 = %d, best t2 = %d\n", i, i-j, tactus[i][j].best_triple1, tactus[i][j].best_triple2);
	}
    }
}

double first_beat_note_func(int p) {

    if(note_pip[p]==1) return log(.6);
    else return log(.4);
}

double note_func(int phase, int p) {

    /* note_pip[p]==1 means there's a note-onset there; note_pip[p]==0 means there isn't */
    /* phase=0 means it's a level 3 beat, phase=1 or 2 means it's a level 2 beat */

    /* If it's the last pip of the piece, we should return 1 - we assume there's not a note there. */

    if(note_pip[p]==1) {
	if(p == last_pip) return log(0.01); /* A note-onset on the very last pip?? Probably never happens */
	else if(phase==0) return log(UW);
	else return log(TW);
    }
    
    else {
	if(p == last_pip) return log(0.99);
	else if(phase==0) return log(1.0-UW);
	else return log(1.0-TW);
    }
}

double make_tables(int timesig, int upper) {
    
    int i, j, k, bestk, bestj, besti;
    double cell_score[14][3];
    double best_score;
    double i_score;
    double another;
    int ph, best_ph;
    double total_mass=0.0, total_phase_mass;
    double adj_previous_score, present_score, adj_present_score;
    double mass;
    double phase_prob;
    int pph;

    /* "another" is the probability of generating another beat */
    another = .95;

    /* In tactus[x][y], x is a pip number, y is a pip interval size for the previous interval;
       the value is the best analysis ending in that interval at that pip. */

    /* A phase zero analysis of tactus (i-j,i) is one in which the i-j beat is assumed to be 
       strong. This was done because it's logical, at the initial tactuses, for phase 0 to
       imply a strong first beat. However, this means that the current beat for phase 0 
       (the one that you're assigning a note score for) is weak. Similarly, in 3/4,
       phase 1 means the i-j beat is beat 2, phase 2 means the i-j beat is beat 3. */

    /* I _think_ if you set upper==1, it will basically generate just levels 1 and 2,
       though it will do this several times. This hasn't really been tested, though. */

    /* We now assume that the first note-onset is at pip 21. The location of the first tactus
       beat might be anywhere between pips 0 and 21, inclusive. 

    /* First generate the initial tables */

    for(i=22; i<44; i++) {
	for(j=9; j<23; j++) {
	    for(ph=0; ph<upper; ph++) {
		
		if(upper==1) phase_prob = 1.0;
		if(upper==2 && ph==0) {
		    if(note_pip[i-j]) phase_prob = DZP;
		    else phase_prob = 1.0-DZP;
		}
		if(upper==2 && ph==1) {
		    if(note_pip[i-j]) phase_prob = 1.0-DZP;
		    else phase_prob = DZP;
		}
		if(upper==3 && ph==0) {
		    if(note_pip[i-j]) phase_prob = TZP;
		    else phase_prob = (1.0-TZP)*(1.0-TTP);
		}
		if(upper==3 && ph==2) {
		    if(note_pip[i-j]) phase_prob = (1.0-TZP)*TTP;
		    else phase_prob = TZP;
		}
		if(upper==3 && ph==1) {
		    if(note_pip[i-j]) phase_prob = (1.0-TZP)*(1.0-TTP);
		    else phase_prob = (1.0-TZP)*TTP;
		}		

		if(timesig==2) {
		    /* This next line (and similar ones below) shouldn't be necessary - the
		       scores and masses for these cells will all be set later. (You could just
		       "continue") */
		    if(i-j > 21) tactus[i][j].score[ph] = -1000.0;
		    else tactus[i][j].score[ph] = log(first_dist[j-9]) + tactus[i][j].best_duple_score + first_beat_note_func(i-j) + note_func((ph+1)%upper,i) + log(phase_prob);
		}
		if(timesig==3) {
		    if(i-j > 21) tactus[i][j].score[ph] = -1000.0;
		    else tactus[i][j].score[ph] = log(first_dist[j-9]) + tactus[i][j].best_triple_score + first_beat_note_func(i-j) + note_func((ph+1)%upper,i) + log(phase_prob);
		}
		/* The line below indicates when something is an initial tactus interval */
		tactus[i][j].bestk[ph]=0;
		
		//if(upper==2 && ph==0) printf("i=%d, j=%d: tactus score = %6.3f\n", i, j, tactus[i][j].score[ph]);

		if(timesig==2) {
		    if(i-j>21) tactus[i][j].mass[ph] = 0.0;
		    else tactus[i][j].mass[ph] = exp((log(first_dist[j-9]) + log(tactus[i][j].duple_mass) + first_beat_note_func(i-j) + note_func((ph+1)%upper,i) + log(phase_prob)) / sec(i));
		}
		if(timesig==3) {
		    if(i-j>21) tactus[i][j].mass[ph] = 0.0;
		    tactus[i][j].mass[ph] = exp((log(first_dist[j-9]) + log(tactus[i][j].triple_mass) + first_beat_note_func(i-j) + note_func((ph+1)%upper,i) + log(phase_prob)) / sec(i));
		}
	    }
	}
    }

    /* Now generate non-initial analyses. We skip all cases where i-j-k<0 (this automatically
       skips all cases where i-j < 9). We go through all the k's twice: once to set the
       scores, the second to choose the best one and to set the masses. We could probably
       do it all on one pass instead. */

    for(i=31; i<=last_pip; i++) {

	for(ph=0; ph<upper; ph++) {

	    for(j=9; j<23; j++) {

		/* We just need the following line so it doesn't go through the second k loop */		
		if(i-j < 22) continue;
		
		tactus[i][j].mass[ph] = 0.0;
		
		for(k=9; k<23; k++) {
		    
		    if(interval_score[j-9][k-9]==0.0) i_score = -1000.0;
		    else i_score = log(interval_score[j-9][k-9]);
		    
		    pph = (upper+ph-1)%upper;    /* pph = previous phase */
		    if(timesig==2) cell_score[k-9][ph] = tactus[i-j][k].score[pph] + tactus[i][j].best_duple_score + i_score + note_func( (ph+1)%upper, i) + log(another);
		    if(timesig==3) cell_score[k-9][ph] = tactus[i-j][k].score[pph] + tactus[i][j].best_triple_score + i_score + note_func( (ph+1)%upper, i) + log(another);
		    //if(i==43 && j==21 && upper==3 && ph==0) printf("i=%d, j=%d, k=%d: previous score = %6.3f; current score = %6.3f\n", i, j, k, tactus[i-j][k].score[pph], cell_score[k-9][ph]);

		}

		best_score = -1000000.0;
		bestk = -1;

		mass = 0.0;
		for(k=9; k<23; k++) {
		    if(cell_score[k-9][ph] > best_score) {
			best_score = cell_score[k-9][ph];
			bestk = k;
		    }

		    //if(i==100) printf("best k = %d\n", bestk);
		    if(bestk == -1) {printf("Error - at i=%d, j=%d, no best k found!\n", i, j); exit(1);}

		    tactus[i][j].bestk[ph] = bestk;
		    tactus[i][j].score[ph] = best_score;

		    if(interval_score[j-9][k-9]==0.0) i_score = -1000.0;
		    else i_score = log(interval_score[j-9][k-9]);
		    pph = ((upper+ph)-1)%upper;

		    if(timesig==2) {

			/* adj_prev_score should be a large negative number; pres should be a probability; adj_pres should be a large negative number; tactus[i][j].mass should be a probability */
			adj_previous_score = log(tactus[i-j][k].mass[pph]) * sec(i-j);
			present_score = exp(i_score) * tactus[i][j].duple_mass * exp(note_func( (ph+1)%upper, i)) * another;
			adj_present_score = log(present_score);
			mass += exp(adj_previous_score + adj_present_score);
			//printf("i=%d, j=%d, sec(i-j)=%6.3f, k=%d: adj_prev = %6.3f, adj_pres = %6.3f, mass += %f\n", i, j, sec(i-j), k, adj_previous_score, adj_present_score, exp((adj_previous_score + adj_present_score)/sec(i)));
			if(v > -1 && adj_previous_score >= 0.0) {printf("Doh!\n"), exit(1); }
			if(v > -1 && present_score > 1.0) {printf("Whoops!\n"), exit(1); }
			if(v > -1 && exp(adj_previous_score + adj_present_score) > 1.0) { printf("Oops!\n"); exit(1); }
		    }
		    if(timesig==3) {
			adj_previous_score = log(tactus[i-j][k].mass[pph]) * sec(i-j);
			present_score = exp(i_score) * tactus[i][j].triple_mass * exp(note_func( (ph+1)%upper, i)) * another;
			adj_present_score = log(present_score);
			mass += exp(adj_previous_score + adj_present_score);
		    }

		}

		tactus[i][j].mass[ph] = pow(mass, (1.0/sec(i)));

		//printf("i=%d, j=%d: mass = %f\n", i, j, tactus[i][j].mass[ph]);

		//if(i==100) printf("Final score for i=%d, j=%d: %6.3f\n", i, j, tactus[i][j].score);
	    }
	}
    }

    /* Now choose the best final j and the best phase. The best phase here represents
     the phase of the final i-j beat (i.e. the penultimate beat), NOT the first beat! */

    best_score = -1000000.0;
    bestj=-1;

    i=last_pip;

    for(ph=0; ph<upper; ph++) {
	total_phase_mass = 0.0;
	for(j=9; j<23; j++) {

	    //printf("j=%d: total mass before final adjustments=%f\n", j, tactus[i][j].mass[ph]);

	    /* Now we add on the P of NOT generating another beat (1.0-another) */
	    tactus[i][j].score[ph] += log(1.0-another);
	    tactus[i][j].mass[ph] *= pow((1.0-another), 1.0/sec(i));

	    if(timesig==2) {
		tactus[i][j].score[ph] += log(LDP);
		tactus[i][j].mass[ph] *= pow(LDP, 1.0/sec(i));
	    }
	    if(timesig==3) {
		tactus[i][j].score[ph] += log(1.0-LDP);
		tactus[i][j].mass[ph] *= pow((1.0-LDP), 1.0/sec(i));
	    }
	    if(upper==2) {
		tactus[i][j].score[ph] += log(UDP);
		tactus[i][j].mass[ph] *= pow(UDP, 1.0/sec(i));
	    }
	    if(upper==3) {
		tactus[i][j].score[ph] += log(1.0-UDP);
		tactus[i][j].mass[ph] *= pow((1.0-UDP), 1.0/sec(i));
	    }
	    //printf("after=%f\n", tactus[i][j].mass);

	    if(tactus[i][j].score[ph] > best_score) {
		best_score = tactus[i][j].score[ph];
		bestj = j;
		best_ph=ph;
	    }
	    
	    total_phase_mass += pow(tactus[i][j].mass[ph], sec(i));
	}
	//printf("Final best score for last i=%d, ph=%d: best j=%d, score=%6.3f\n", i, ph, bestj, tactus[i][bestj].score);
	if(v>1) printf("    log phase mass for ph=%d: %6.3f\n", ph, log(total_phase_mass));
	total_mass += total_phase_mass;
    }

    grand_total_mass += total_mass;

    if(bestj == -1) {printf("Error - no best final j found!\n"); exit(1);}

    if(v>1) printf("  Total mass = %6.40f, log=%6.3f\n", total_mass, log(total_mass));

    //printf("Best phase = %d\n", best_ph);

    best_final_j[timesig-2][upper-2] = bestj;
    best_phase[timesig-2][upper-2] = best_ph;

    if(v>1) printf("  log P of best analysis = %6.3f\n", best_score);
    return best_score;

}

void trace_back(int timesig, int upper) {

    int i, j, newj, p, ph, previous_tactus, count, phase, num_tactus_beats, first_tactus;
    int final_phase[10000];

    /* Now we set the final beat strength and final phase of each
        tactus beat. The phase of a beat indicates its actual phase,
        i.e. phase 0 means it's a beat at level 3, phase 1 means it's
        one after a beat at level 3, etc.. */

    for(p=0; p<=last_pip; p++) {
	final[p]=0;
	final_phase[p]=-1;
    }

    //best_phase[0][0] = 1;

    i=last_pip;
    j=best_final_j[timesig-2][upper-2];
    ph = phase = best_phase[timesig-2][upper-2];

    final_phase[i] = (phase+1) % upper;
    final_phase[i-j] = phase;

    while(1) {
	newj = tactus[i][j].bestk[ph];
	i = i-j;
	j = newj;
	ph = (upper+ph-1)%upper;
	final_phase[i-j] = ph;
	//printf("%d ", final_phase[i-j]);
	if(tactus[i][j].bestk[ph] == 0) break;
    }
    //printf("\n");

    for(p=0; p<=last_pip; p++) {
	if(final_phase[p]==-1) continue;
	else if(final_phase[p]==0) final[p]=3;
	else final[p]=2;
    }

    /* Now we define the lower level beats */

    previous_tactus = -1;
    for(p=0; p<=last_pip; p++) {
	if(final[p]>0) {
	    if(previous_tactus == -1) {
		first_tactus = p;
		previous_tactus = p;
		continue;
	    }
	    if(timesig==2) final[tactus[p][p-previous_tactus].best_duple]=1;
	    if(timesig==3) {
		final[tactus[p][p-previous_tactus].best_triple1]=1;
		final[tactus[p][p-previous_tactus].best_triple2]=1;
	    }	    
	    previous_tactus = p;
	}
    }

    if(v>1) printf("Best initial phase = %d\n", final_phase[0]);
    for(p=0, num_tactus_beats=0; p<=last_pip; p++) if(final[p]>1) num_tactus_beats++;
    if(v>1) printf("%d tactus beats, average length = %d\n", num_tactus_beats, 
	   (last_pip * 50) / num_tactus_beats);

    /* Terse display, just beat positions and phases */
    if(v>1) {
	for(p=0; p<=last_pip; p++) {
	    if(final[p] > 1) {
		printf("%d", p);
		printf("[%d] ", final_phase[p]);
	    }
	    else if(final[p]==1) printf("(%d) ", p);
	}
	printf("\n");
    }

    /* Verbose display, showing scores etc. */

    /*
    previous_tactus = first_tactus;
    for(p=first_tactus+1; p<=last_pip; p++) {
	if(final[p] > 1) {
	    j = p-previous_tactus;
	    ph = (upper+final_phase[p]-1)%upper;
	    printf("i=%d, j=%d, ph = %d, level=%d, best k = %d, score=%6.3f\n", p, j, ph, final[p], tactus[p][j].bestk[ph], tactus[p][j].score[ph]);
	    previous_tactus = p;
	}
	}  */

}

int print_note(int actually_print, int pip) {

    int n, pc;
    for(n=0; n<numnotes; n++) {
	if(note[n].ontime / 50 == pip) {
	    pc = note[n].pitch % 12;
	    if(actually_print) {
		if(pc == 0) printf("C ");
		else if(pc == 1) printf("C# ");
		else if(pc == 2) printf("D ");
		else if(pc == 3) printf("Eb ");
		else if(pc == 4) printf("E ");
		else if(pc == 5) printf("F ");
		else if(pc == 6) printf("F# ");
		else if(pc == 7) printf("G ");
		else if(pc == 8) printf("Ab ");
		else if(pc == 9) printf("A ");
		else if(pc == 10) printf("Bb ");
		else if(pc == 11) printf("B ");
		else printf("? ");
	    }
	    if(pc == 1 || pc == 3 || pc == 6 || pc == 8 || pc == 10) return 2; 
	    else return 1;
	}
    }
    /* We should never get here */
    return 1;
}

void print_graphic_output() {

    int np, np2, level, first_onset=-1;

    np=0;

    for(np=0; np<=last_pip; np++) {
	if(note_pip[np]) {
	    first_onset = np;
	    break;
	}
    }
    //printf("First onset = %d\n", first_onset);
    last_pip -= first_onset;
    
    if(first_onset != -1) {
	for(np=0; np<=last_pip; np++) {
	    final[np] = final[np+first_onset];
	    note_pip[np] = note_pip[np+first_onset];
	}
    }

    np=0;
    while(1) {

        for(level = 3; level>=1; level--) {

	    printf("\n     ");
	    for(np2=0; np2<40; np2++) {
	      if(np+np2 > last_pip) break;
	      if(np2%10 == 0) printf("  ");
	      if(final[np+np2] >= level) printf("x ");
	      else printf("  ");
	      if(note_pip[np+np2]==1) if(print_note(0, (np+np2+first_onset))==2) printf(" ");
	    }
	}

	printf("\n%4d ", np);
	for(np2=0; np2<40; np2++) {
	    if(np2%10 == 0) printf("| ");
	    if(np+np2 > last_pip) break;
	    if(note_pip[np+np2]==1) print_note(1, (np+np2+first_onset));
	    else printf("o ");
	}
	if(np+np2 > last_pip) break;
	np += 40;
    } 
    printf("\n");
}

void print_notebeat_list() {

    int p, n, adj_ontime, adj_offtime, beat_time;
    for(p=0; p<=last_pip; p++) {
	if(final[p] > 0) {
	    beat_time = (p*50) - note_offset;
	    if(beat_time < 0) continue;
	    printf("Beat %d %d\n", beat_time, final[p]);
	}
    }

    for(n=0; n<numnotes; n++) {
	adj_ontime = ((note[n].ontime / 50) * 50) - note_offset;
	adj_offtime =  ((note[n].offtime / 50) * 50) - note_offset;
	printf("Note %d %d %d\n", adj_ontime, adj_offtime, note[n].pitch);
    }
}


