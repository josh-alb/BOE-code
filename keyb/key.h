/* See main.c for general info */

FILE *in_file;
FILE *out_file;

int display_command;
int numnotes, numbeats, numchords, num_sbeats;
int total_duration;
int final_timepoint;
int firstbeat;
int harmonic_input;
extern double default_profile_value;
double major_profile[12]; /* In line of fifths order */
double minor_profile[12];
extern int segment_beat_level;
extern int beat_printout_level;
extern float change_penalty;
extern int verbosity;
extern int romnums;
extern int romnum_type;
extern int running;
extern int npc_or_tpc_profile;
extern int scoring_mode;
double seglength;
typedef struct note_struct {
  int ontime;
  int offtime;
  int duration;
  int pitch;
  int tpc;
} blah;

struct note_struct note[10000];

typedef struct chord_struct {
  int ontime;
  int offtime;
  int duration;
  int root;
  int key;
  int function;
  int romnum;
  int mode;       /* major =1, minor =2, unspecified = 0 */
  int ext;        /* seventh = 1, no seventh = 0 */
  int inversion;  /* 0 = root, 1 = first, 2 = second, 3 = third */
  int fifth;      /* 1 = perfect fifth, 2 = diminished, 3 = unspecified */
  int beatlevel;
} blar;

struct chord_struct ichord[5000];      /* input chords */
struct chord_struct chord[2000];       
struct {
  int time;
  int level;
} beat[5000];

struct {
  int time;
} sbeat[1000];

int pc_dur[12];
int pc_tally[500];

extern char letter[7];

double key_profile[56][28];

struct {
    int start;
    int end;
    struct note_struct snote[100];
    int numnotes;        /* number of notes in the segment */
    double average_dur;         /* average input vector value (needed for K-S algorithm) */
} segment[500];        /* An array storing the notes in each segment. */
int segtotal;              /* total number of segments - 1 */

int seg_prof[500][28];
double key_score[500][56];
double total_prob[500];

double first_seg_analysis[56];

double analysis[500][56][56];
int best[500][56];
extern int seg;
int final[500];
int provisional[500][500]; /* used for printing out provisional analyses at each step */

double keyfit_score[500];

int main(int, char **);
    
void bad_param(char *);
void read_parameter_file(char *, int);
void bad_input(int);

void print_keyname(int);
void create_chords();
void create_segments();
void fill_segments();
void count_segment_notes();
void prepare_profiles();
void generate_tpc_profiles();
void generate_npc_profiles();
void match_profiles();
void make_first_table();
void choose_best_i();
void make_tables();
void best_key_analysis();

void generate_chord_info();
void merge_functions();
void chords_to_romnums_kp();
void print_romnums_kp();
void chords_to_romnums_as();
void print_romnums_as();
void display_table();
void display_running();

