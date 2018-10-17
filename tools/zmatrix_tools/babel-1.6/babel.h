#ifndef __BABEL_BABEL_H__
#define __BABEL_BABEL_H__

#ifdef __cplusplus
extern "C" {
#endif

/* addh.c */
extern void add_hydrogens(ums_type *mol);
extern void place_hydrogens1(ums_type *mol, int old_count, int num_H_to_add);
extern void add_methyl_hydrogen(ums_type *mol, int c_num, int h_num, double b_length);
extern void add_sp_hydrogen(ums_type *mol, int c_num, int h_num, double b_length);
extern void add_sp2_hydrogen(ums_type *mol, int c_num, int h_num, double b_length);
extern void add_vinyl_hydrogens(ums_type *mol, int c_num, int h_num, double b_length);
extern void add_sp3_N_hydrogen(ums_type *mol, int n_num, int h_num, double b_length);
extern void add_methylene_hydrogens(ums_type *mol, int c_num, int h_num, double b_length);
extern void add_tertiary_hydrogen(ums_type *mol, int c_num, int h_num, double b_length);
extern int type_added_hydrogen(ums_type *mol, int c1, int h1);
extern int count_missing_hydrogens(ums_type *mol);
extern void add_2d_hydrogens(ums_type *mol);
/* addh2.c */
extern int count_missing_bo_hydrogens(ums_type *mol);
extern void place_hydrogens2(ums_type *mol, int old_count, int num_H_to_add);
extern int count_attached_bonds(ums_type *mol, int atm);
/* aromatic.c */
extern void find_aromatic_atoms(ums_type *mol);
extern int count_arom_atm_bonds(ums_type *mol, int atm, ring_info *info, int which_ring);
extern void setup_ring_info(ums_type *mol, ring_struct *rings, ring_info *info);
extern void cleanup_ring_info(ums_type *mol, ring_struct *rings, ring_info *info);
extern void print_ring_info(ums_type *mol, ring_info *info);
extern void find_aromatic_rings2(ums_type *mol, ring_struct *rings, ring_info *info);
extern void find_aromatic_rings(ums_type *mol, ring_struct *rings, ring_info *info);
extern int bond_is_aromatic(int a, int b, ring_info *info);
extern int in_same_ring(int a, int b, ring_info *info, int which_ring);
extern int count_ortho_substituents(ums_type *mol, int start, int end);
extern int count_bonded_hydrogens(ums_type *mol, int atm);
/* assbnd.c */
extern int assign_radii(ums_type *mol);
extern void fast_radii(ums_type *mol);
extern int assign_hybrid_radii(ums_type *mol);
extern int assign_bonds(ums_type *mol);
extern int assign_pdb_bonds(ums_type *mol);
extern void estimate_bond_order(ums_type *mol);
extern void check_bonds(ums_type *mol);
extern int find_close_atoms(ums_type *mol, int num);
extern int zsort_atoms(temp_atom_rec *a, temp_atom_rec *b);
extern int sort_connections(connect_type *a, connect_type *b);
extern double distance(coord_type first, coord_type second);
extern int is_element(char *ele);
extern int get_atomic_number(char *sym);
/* asstypes.c */
extern void check_atomic_numbers(ums_type *mol);
extern int assign_types(ums_type *mol);
extern void tag_organics(ums_type *mol);
extern int phase1(ums_type *mol);
extern int phase2(ums_type *mol);
extern int phase3(ums_type *mol);
extern int phase4(ums_type *mol);
extern int phase5(ums_type *mol);
extern int phase6(ums_type *mol);
extern int type_hydrogens(ums_type *mol);
extern int valence_four(ums_type *mol);
extern int valence_three(ums_type *mol);
extern int valence_two(ums_type *mol);
extern int valence_one(ums_type *mol);
extern double bond_angle(coord_type a, coord_type b, coord_type c);
extern int count_heavy_atoms(ums_type *mol, int atom_number);
extern int count_free_ox(ums_type *mol, int atom_number);
extern void fix_carboxylates(ums_type *mol);
extern void check_for_amides(ums_type *mol);
/* block.c */
extern int block_alloc(block_ptr *handle, const char *types, ...);
extern int block_calloc(block_ptr *handle, const char *types, ...);
extern void block_free(block_ptr *handle);
/* bndord.c */
extern void assign_bond_order(ums_type *mol);
extern void assign_bond_order1(ums_type *mol);
extern int is_carboxyl(ums_type *mol, int the_bond);
extern int check_for_overflow(ums_type *mol);
extern int check_for_conjugation(ums_type *mol);
extern int atom_in_common(ums_type *mol, int bond1, int bond2);
extern int assign_bond_code(char *input_type);
/* bo.c */
extern void assign_bond_order2(ums_type *mol);
extern void set_bo(ums_type *mol);
extern int connect_the_dots(ums_type *mol, int atm, int start, int *dots, bnd_stack *stk);
extern int get_bond_number(ums_type *mol, int start, int end);
extern void tag_ring_atoms(ums_type *mol, path *ring_set, int ring_count);
extern double get_bond_ratio(ums_type *mol, int a1, int a2);
extern void estimate_bond_order2(ums_type *mol);
extern void process_5_ring(ums_type *mol, ring_struct *rings, int num);
extern void print_bo(ums_type *mol);
extern void reset_bonds(ums_type *mol, int *temp_bo);
extern void save_bond_orders(ums_type *mol, int *temp_bo);
extern int check_bond_order(ums_type *mol);
extern double torsion_angle(coord_type a, coord_type b, coord_type c, coord_type d);
extern void dearomatize(ums_type *mol);
extern int has_aromatic_bonds(ums_type *mol);
/* buildct.c */
extern int build_connection_table(ums_type *mol);
extern int member(int first, int second, connect_type *connections, int length);
/* combine.c */
extern int borrow_types(ums_type *mol1, ums_type *mol2);
/* convert.c */
extern void new_control(ums_type *mol);
extern void init_babel_control(ums_type *mol);
extern void process_command_line_args(int argc, char *argv[], ums_type *mol);
extern void translate_input_code(char *code, ums_type *mol);
extern void translate_output_code(char *code, ums_type *mol);
extern ums_type *do_inputs(ums_type *mol);
extern ums_type *do_outputs(FILE *outfile, ums_type *mol);
extern void generate_outfile_name(ums_type *mol, int count);
extern void set_limits(ums_type *mol);
extern int want_this_file(ums_type *mol, int counter, int end);
extern void usage(void);
extern void write_bbldef(void);
extern void show_inputs(void);
extern void show_outputs(void);
extern int output_all_formats(FILE *fp, ums_type *mol);
extern void write_blurb(void);
/* delatms.c */
extern ums_type *delete_atoms(ums_type *mol, char *del_str);
extern int tag_atoms(ums_type *mol, char *del_str);
extern ums_type *build_new_ums(ums_type *mol, int heavy_count);
extern void dissect_connection_table(ums_type *mol);
extern void translate_del_str(char *inp_type, char *del_str);
/* delh2o.c */
extern ums_type *delete_water(ums_type *mol);
extern int tag_waters(ums_type *mol);
extern int is_water(ums_type *mol, int i);
/* filesrch.c */
extern long ffsearch(FILE *fp, const char *pattern, const size_t size, int N);
extern long rfsearch(FILE *fp, const char *pattern, const size_t size, int N);
/* fileutil.c */
extern FILE *open_w_env(char *f_name, char *env_var);
extern FILE *open_read(char *file_name);
extern FILE *open_write(char *file_name);
extern void close_file(char *file_name, FILE *file1);
/* gastchg.c */
extern void calc_gasteiger_charges(ums_type *mol);
extern void setup_sigma_params(gast_param *par);
extern void lookup_sigma_params(ums_type *mol, gast_param *master_par, gast_param *atm_par);
extern void calc_sigma_charges(ums_type *mol);
extern void print_gasteiger_params(ums_type *mol, gast_param *atm_par);
extern int is_carboxylate(ums_type *mol, int atm);
/* htoend.c */
extern void shift_h_to_end(ums_type *mol);
/* int2cart.c */
extern int old_int_to_cart(ums_type *mol);
/* intcart.c */
extern int int_to_cart(ums_type *mol);
/* menus.c */
extern void babel_menus(ums_type *mol);
extern void read_menu(ums_type *mol);
extern void write_menu(ums_type *mol);
extern int continuation_menu(ums_type *mol);
extern int get_choice(int min, int max);
/* miniums.c */
extern void ums_to_mini_ums(mini_ums *mini, ums_type *mol);
extern void copy_coordinates(mini_ums *mini, ums_type *mol);
extern void copy_mini(mini_ums *mini1, mini_ums *mini2);
extern void initialize_mini(mini_ums *mini);
extern void make_mini_enantiomer(mini_ums *ena, mini_ums *mini);
extern void print_mini(mini_ums *mini);
extern void add_mini(mini_ums *mini1, mini_ums *mini2, mini_ums *new_mini);
extern void mini_to_set(mini_ums *mini, set_type *the_set);
extern int non_zero(coord_type *pt);
extern void adjust_mini_vector(mini_ums *mini, int core_atom, int vect_atom, double new_dist);
extern void release_mini(mini_ums *mini);
extern void zero_mini(mini_ums *mini);
extern void read_mini(FILE *file, mini_ums *mini);
extern long int write_mini(FILE *file, mini_ums *mini);
/* molwt.c */
void get_atomic_weights(ums_type *mol);
double calc_molecular_weight(ums_type *mol);
double atnum2mass(int atnum);
extern double atnum2mass(int atnum);
/* nodummy.c */
extern int dummy_check(ums_type *mol);
/* orient.c */
extern void set_std_orientation(ums_type *mol);
extern void get_orientation_matrix(ums_type *mol, int x1, int x2, int y1, int y2, matrix_3x3 *m);
extern void ums_plus_vector(ums_type *mol, vect_type *v);
extern void ums_dot_matrix(ums_type *mol, matrix_3x3 *m);
extern void orient_ums(ums_type *mol, int x1, int x2, int y1, int y2);
extern void coord_to_vector(vect_type *v, coord_type *c);
extern void vector_to_coord(coord_type *c, vect_type *v);
extern void center_at_atom(ums_type *mol, int atom);
extern void center_at_atoms(ums_type *mol, int atom[], int count);
extern void center_at_origin(ums_type *mol, vect_type *v);
extern void AlignMol(ums_type *mol);
void find_mol_center(ums_type *mol, vect_type *center, int use_wts);

/* precip.c */
extern ums_type *precipitate(ums_type *mol);
extern int tag_salts(ums_type *mol);
extern void copy_ums(ums_type *dest, ums_type *src);
/* printbad.c */
extern void print_bad_connections(ums_type *mol, int atom);
/* progress.c */
extern void ShowProgress(int i, char *the_string);
extern void UpdateProgress(void);
/* rdalch.c */
extern int read_alchemy(FILE *file1, ums_type *mol);
extern int translate_alchemy_bond_order(char *bo_string);
/* rdampout.c */
extern int read_ampac_output(FILE *file1, ums_type *mol);
/* rdbalst.c */
extern int read_bs(FILE *file1, ums_type *mol);
/* rdbgf.c */
extern int read_bgf(FILE *file1, ums_type *mol);
/* rdboogie.c */
extern int read_boogie(FILE *file1, ums_type *mol);
/* rdc3d.c */
extern int read_mmads(FILE *file1, ums_type *mol);
extern int read_chem3d1(FILE *file1, ums_type *mol);
extern int read_chem3d2(FILE *file1, ums_type *mol);
extern int read_chem3d(FILE *file1, ums_type *mol, char *keywords, char *type_key);
extern void xlate_c3d_labels(ums_type *mol, int *label);
/* rdcacao.c */
extern int read_caccrt(FILE *file1, ums_type *mol);
/* rdcadpac.c */
extern int read_cadpac(FILE *file1, ums_type *mol);
/* rdcharmm.c */
extern int read_charmm(FILE *file1, ums_type *mol);
/* rdcsd.c */
extern int read_csd(FILE *file1, ums_type *mol);
/* rddock.c */
extern int read_dock_database(FILE *file1, ums_type *mol);
/* rddpdb.c */
extern int read_dock_pdb(FILE *file1, ums_type *mol);
/* rdelmnts.c */
extern int read_element_file(void);
extern void clean_up_elements(void);
extern int write_elements(element_type elements[]);
extern int get_max_bonds(int atomic_number);
extern void atomic_number_to_name(int i, char *name);
extern void add_element_types(ums_type *mol);
/* rdfdat.c */
extern int read_fdat(FILE *file1, ums_type *mol);
extern void my_strncpy(char *str1, char *str2, int len);
extern double my_atof(char *the_str);
extern double my_atoi(char *the_str);
extern int check_refcode(char *keywords, char *csd_line, ums_type *mol);
/* rdfeat.c */
extern int read_feat(FILE *file1, ums_type *mol);
/* rdfract.c */
extern int read_csd_fractional(FILE *file1, ums_type *mol);
extern int read_fform_fract(FILE *file1, ums_type *mol);
extern void fill_orth_matrix(fract_type *f, matrix_3x3 *m);
extern void fract_to_cart(coord_type *p, matrix_3x3 *m);
/* rdgamout.c */
extern int read_gamess_output(FILE *file1, ums_type *mol);
extern void bohr_to_angstroms(ums_type *mol);
/* rdgauout.c */
extern int read_gau_out(FILE *file1, ums_type *mol);
extern void read_gaussian_cartesian_atoms(FILE *file1, ums_type *mol);
extern void read_gaussian_internal_atoms(FILE *file1, ums_type *mol);
/* rdgzmat.c */
extern int read_gau_zmatrix(FILE *file1, ums_type *mol);
extern int xlate_label(ums_type *mol, char *sym, int max);
extern double xlate_symbol(zsymbol *the_syms, char *sym, int max);
/* rdhin.c */
extern int read_hin(FILE *file1, ums_type *mol);
extern int translate_hin_bond_order(char c);
/* rdinsite.c */
extern int read_biosym_car(FILE *file1, ums_type *mol);
/* rdint.c */
extern int read_mopint(FILE *file1, ums_type *mol);
/* rdirc.c */
extern int read_gaussian_94(FILE *fp, ums_type *mol);
/* rdisis.c */
extern int read_isis(FILE *file1, ums_type *mol);
/* rdm3d.c */
extern int read_m3d(FILE *file1, ums_type *mol);
extern int translate_m3d_bond_order(char *bo_string);
/* rdmacmod.c */
extern int read_macromodel(FILE *file1, ums_type *mol);
extern void figure_valence(ums_type *mol);
extern int assign_mmd_bond_order(ums_type *mol);
/* rdmacmol.c */
extern int read_mcmol(FILE *file1, ums_type *mol);
extern char *clean_comments(char *mcmol_str);
/* rdmdl.c */
extern int read_mdl(FILE *file1, ums_type *mol);
extern double get_scale_factor(ums_type *mol);
extern void scale_for_ChemWindow(ums_type *mol);
/* rdmicro.c */
extern int read_micro(FILE *file1, ums_type *mol);
extern char *strip_front_num(char *the_str);
/* rdmm2.c */
extern int read_mm2(FILE *file1, ums_type *mol);
/* rdmm2in.c */
extern int read_mm2_input(FILE *file1, ums_type *mol);
/* rdmm3.c */
extern int read_mm3(FILE *file1, ums_type *mol);
/* rdmolen.c */
extern int read_molin(FILE *file1, ums_type *mol);
extern void get_cell_params(fract_type *f);
/* rdmopac.c */
extern int read_mopac_output(FILE *file1, ums_type *mol);
/* rdmopcrt.c */
extern int read_mop_cart(FILE *file1, ums_type *mol);
/* rdpcmod.c */
extern int read_pcmodel(FILE *file1, ums_type *mol);
extern void get_pcmod_bonds(char *the_line, ums_type *mol, int i);
/* rdpdb.c */
extern int read_pdb(FILE *file1, ums_type *mol);
extern void fix_A_type(char *type, char *id, char *res_type);
extern void process_connect_records(FILE *file1, ums_type *mol);
extern void add_bond(ums_type *mol, int i, int j);
extern int parse_atom_record(char *pdb_line, int the_atom, ums_type *mol);
extern int parse_conect_record(char *pdb_line, int *serial_number,
  int *conn1, int *conn2, int *conn3, int *conn4);
extern int stringToType(int record_type, int the_atom, ums_type *mol);
/* rdprep.c */
extern int read_amber_prep(FILE *file1, ums_type *mol);
/* rdpsgout.c */
extern int read_psgvb_output(FILE *file1, ums_type *mol);
/* rdpsgvin.c */
extern int read_psgvb_input(FILE *file1, ums_type *mol);
extern double xlate_ps_symbol(zsymbol *the_syms, char *sym, int max);
/* rdquanta.c */
extern int read_quanta(FILE *file1, ums_type *mol);
extern int read_quanta_types(char3 *quanta_types);
/* rdschak.c */
extern int read_schakal(FILE *file1, ums_type *mol);
/* rdshelx.c */
extern int read_shelx(FILE *file1, ums_type *mol);
extern int count_shelx_atoms(char *the_line);
extern void check_shelx_coords(coord_type *p);
extern int is_good_shelx_line(char *the_line);
/* rdsmiles.c */
extern int read_smiles(FILE *fp, ums_type *mol);
extern void ct_to_ums(smilescontab_t *ct, ums_type *mol);
/* rdspart.c */
extern int read_spartan(FILE *file1, ums_type *mol);
/* rdspmm.c */
extern int read_spartan_mol_mech(FILE *file1, ums_type *mol);
/* rdspsemi.c */
extern int read_spartan_semiempirical(FILE *file1, ums_type *mol);
/* rdsybmol.c */
extern int read_sybyl_mol(FILE *file1, ums_type *mol);
/* rdsybyl.c */
extern int read_sybyl(FILE *file1, ums_type *mol);
/* rdtypes.c */
extern void read_types_table(void);
extern void clean_up_types_table(void);
extern void write_type_table(void);
extern int locate_input_type(char *format);
extern int get_input_type(int atom_num, int col_num, char *input, char *std_type, enum type_err error);
extern int get_output_type(int at_num, char *format, char *input, char *out_type, enum type_err error);
extern int get_std_type(char *format, char *input, char *std_type);
extern int xlate_std_type(char *format, char *std_type, char *output);
/* rdunichm.c */
extern int read_unichem(FILE *file1, ums_type *mol);
/* rdwiz.c */
extern int read_wizard(FILE *file1, ums_type *mol);
extern void OrderConnections(ums_type *mol);
/* rdxed.c */
extern int read_xed(FILE *file1, ums_type *mol);
extern void xed_add_connection(ums_type *mol, int start, int end, int bnd_num);
/* rdxyz.c */
extern int read_xyz(FILE *file1, ums_type *mol);
/* renum.c */
extern ums_type *renumber(ums_type *mol);
extern void find_dist_from_origin(ums_type *mol);
extern void sort_by_dist_to_origin(ums_type *mol);
extern int sort_by_dist(temp_atom_rec *a, temp_atom_rec *b);
/* report.c */
extern int print_report_file(FILE *file1, ums_type *mol);
extern void get_element_type(ums_type *mol, int i, char *type);
extern void distance_matrix(ums_type *mol, FILE *file1);
extern void print_torsions(ums_type *mol, FILE *file1);
extern void print_angles(ums_type *mol, FILE *file1);
extern void sort_values(int *a, int *b);
extern int compare_torsion(torsion_rec *t1, torsion_rec *t2);
extern int compare_angle(angle_rec *t1, angle_rec *t2);
/* rings.c */
extern int find_rings(ums_type *mol, ring_struct *rings);
extern void tree_to_sets(ums_type *mol, spanning_tree *the_tree, set_type *the_set[]);
extern void find_common_atom(ums_type *mol, set_type *set1[], set_type *set2[], set_type *set3);
extern void print_tree_set(ums_type *mol, set_type *the_set[]);
extern int sort_rca(connect_type *a, connect_type *b);
extern void build_spanning_tree(ums_type *mol, int root, spanning_tree *tree);
extern void build_restricted_tree(ums_type *mol, int root, int other, spanning_tree *tree);
extern void print_spanning_tree(ums_type *mol, spanning_tree *tree);
extern int find_closure_bonds(ums_type *mol, spanning_tree *tree, connect_type *rca);
extern void path_to_root(spanning_tree *tree, int atom, path *the_path);
extern void init_path(int atoms, path *the_path);
extern void print_path(path *the_path);
extern void paths_to_ring(path path1, path path2, path *the_ring, int atoms);
extern void show_rings(path *ring_set, int ring_count);
extern void sort_rings(path *ring_set, int ring_count);
extern int comp_rings(path *p1, path *p2);
extern void find_bogus_rings(path *ring_set, int ring_count, int atoms);
extern void find_bogus_rings2(ums_type *mol, path *ring_set, int ring_count, int frj);
extern void show_ring(path ring_path);
extern int is_good_ring(path new_ring, path *ring_list, int ring_count);
extern int build_common_array(set_type *common_set, spanning_tree *the_tree, spanning_tree *common_array);
extern int sort_common(spanning_tree *a, spanning_tree *b);
extern void make_ring_ums(ums_type *mol, ums_type *new_mol, path rng);
extern void save_good_rings(ring_struct *good, ring_struct *bad, int count, int dupe_ck);
/* ringutil.c */
extern int find_SSSR(ums_type *mol, ring_struct *rings);
extern ums_type *dissect_ums(ums_type *mol);
extern void add_rings_to_list(ring_struct *rings, ring_struct *tmp, ums_type *mol);
extern ums_type *set_to_ums(ums_type *mol, set_type *set);
extern int find_last_atom(path *the_path);
extern void cleanup_rings(ring_struct *rings);
extern void preserve_rings(ring_struct *good, ring_struct *bad, int count);
/* sets.c */
extern void setclear(set_type *seta);
extern void setcopy(set_type *seta, set_type *setb);
extern void setprint(set_type *seta, char *string);
extern void setand(set_type *seta, set_type *setb, set_type *setc);
extern void setor(set_type *seta, set_type *setb, set_type *setc);
extern void setorxor(set_type *seta, set_type *setb, set_type *setc);
extern void setxor(set_type *seta, set_type *setb, set_type *setc);
extern void setnot(set_type *seta, set_type *setb);
extern int setcmp(set_type *seta, set_type *setb);
extern int setissubset(set_type *seta, set_type *setb);
extern int setcount(set_type *seta);
extern int nextbit(set_type *seta, int last);
extern int NextBit(set_type *seta, int last);
extern set_type *init_set_minbits(int minbits);
extern set_type *init_set_setlen(int setlen);
extern set_type *realloc_set_setlen(set_type *set, int setlen);
extern void free_set(set_type *set);
/* smilesto.c */
extern smilescontab_t *smilestocontab(char *smiles);
extern int addatomtocontab(smilescontab_t **frag, char *atname, int *lastatom, smilesbond_t bondtype);
extern int closecontabring(smilescontab_t *frag, int atom1, int atom2);
extern int aromaticsmilessym(char *sym);
extern int nextfreebondto(smilesatom_t *aptr);
extern int strutils_noccurrences(char *buffer, char c);
/* spline.c */
extern int do_spline(ums_type *mol);
extern void get_vectors(ums_type *start, ums_type *end, vect_type *vect[], int increment);
extern void add_step(ums_type *mol, vect_type *vect[]);
/* strngutl.c */
extern char *rjust(char *str);
extern char *rtrim(char *str);
extern char *ltrim(char *str);
extern char *fill_space(char *dest, char *in, int gesamtl);
extern char *strrev(char *str);
extern int substring_count(char *needle, char *haystack);
extern char *replace_string(char *Str, char *OldStr, char *NewStr);
extern void replace_all_string(char *Str, char *OldStr, char *NewStr);
extern char *ljust(char *str);
/* tokenst.c */
extern int found_token(char *strng, char *delimstr, int zindex);
extern int count_tokens(char *tokens, char *delimstr);
extern char *gettoken(char *tokens, char *delimstr, int tokenindex);
extern void get_token(char *target, char *source, char *delimstr, int tokenindex);
extern char *trim_spaces(char *string);
/* tosmiles.c */
extern char *contabtosmiles(smilescontab_t *frag);
extern int sortranks(int natoms, int *value, int *rank);
extern int resolveties(int natoms, int *oldrank, int *primesum, int *rank);
extern char *generatesmilesstring(smilescontab_t *frag, int *rank);
extern int addsmilesatom(char *string, char *symbol, int aromatic);
extern int addsmilesring(char *string, int ringnum);
extern int lowestunvisitedatom(int natoms, int *rank, int *visited);
extern int findunvisitedring(smilescontab_t *frag, int start, int branch, int *visited);
extern int dotheymeet(smilescontab_t *frag, int start, int next, int target, int *visited);
extern int nunvisitedbondsto(smilesatom_t *aptr, int *visited);
extern int nbondsto(smilesatom_t *aptr);
extern smilesbond_t getbondorder(smilesatom_t *aptr, int atom);
extern int hasaromaticbonds(smilesatom_t *aptr);
/* tree.c */
extern ums_type *renum_for_zmat(ums_type *mol, int base);
extern void build_z_tree(ums_type *mol, int root, z_tree *tree);
extern void push_hydrogens_to_end(ums_type *mol);
extern void continuity_check(ums_type *mol);
extern void find_z_kids(ums_type *mol, z_tree *tree);
extern void print_z_tree(ums_type *mol, z_tree *tree);
extern void dfs_z_tree(int x, ums_type *mol, z_tree *tree);
extern ums_type *renumber_ums(ums_type *mol, int heavy_count);
/* typbybo.c */
extern void assign_type_by_bo(ums_type *mol);
extern void type_attached_oxygens(ums_type *mol, int atm);
/* umslist.c */
extern ums_type *add_ums_to_list(ums_type *new_node, ums_type *list);
extern int count_ums_list(ums_type *list);
extern void cleanup_ums_lst(ums_type *base);
extern ums_type *ums_list_from_file(FILE *fp, int (*reader)());
extern void ums_list_to_file(FILE *fp, ums_type *list, int (*writer)());
extern void show_ums_list(ums_type *list);
/* utils.c */
extern void uppercase(char *str);
extern void lowercase(char *str);
extern void babel_init(void);
extern void babel_cleanup(void);
extern void strip_return(char *the_str);
extern int print_ums(ums_type *mol);
extern void show_warning( const char *format, ... );
extern void fatal_error( const char *format, ... );
extern int print_internal(ums_type *mol);
extern int initialize_ums(ums_type **mol);
extern void zero_out_ums(ums_type *mol, int start);
extern int reinitialize_ums(ums_type **mol);
extern int initialize_internal(ums_type **mol);
extern int initialize_fractional(ums_type **mol);
extern int initialize_residues(ums_type **mol);
extern int release_ums(ums_type *mol);
extern void clean_atom_type(char id[]);
extern double torsion(coord_type a, coord_type b, coord_type c, coord_type d);
extern void free_line(char *the_line);
extern int check_for_eof(FILE *file1);
extern void read_to_eof(FILE *file1);
extern void make_output_name(char *file_name, char *out_name, int counter);
extern void usage_multi_struct(void);
extern void toss(FILE *file1, int num);
extern int is_blank_line(char *the_line);
extern char *new_extension(char *filename, char *extension);
extern int is_one_three(ums_type *mol, int a, int b);
extern double out_of_plane(ums_type *mol, int atm);
/* vectors.c */
extern void point_to_vect(coord_type pt, vect_type *v);
extern void pts_2_vect(ums_type *mol, vect_type *vect, int pt1, int pt2);
extern double vect_ang(vect_type *vect1, vect_type *vect2);
extern double dot(vect_type *vect1, vect_type *vect2);
extern void cross_prod(vect_type *vect1, vect_type *vect2, vect_type *normal);
extern double magnitude(vect_type *vect1);
extern void normalize_vect(vect_type *v1);
extern void vect_sum(vect_type *vect1, vect_type *vect2, vect_type *vect_sm);
extern void vect_diff(vect_type *vect1, vect_type *vect2, vect_type *vect_sm);
extern void scal_x_vect(vect_type *vect1, float scalar);
extern coord_type point_plus_vector(coord_type *p1, vect_type *v1);
extern coord_type point_times_vector(coord_type *p1, vect_type *v1);
extern double determinant_3x3(matrix_3x3 *m);
extern void invert_vector(vect_type *v);
extern void invert_3x3(matrix_3x3 *m);
extern void dump_3x3(matrix_3x3 *m);
extern void mat_3x3_dot_vect(matrix_3x3 *m, vect_type *v);
/* wralch.c */
extern int write_alchemy(FILE *file1, ums_type *mol);
extern void strip_extension(char *file_name, char *new_name);
/* wrbalst.c */
extern int write_bs(FILE *file1, ums_type *mol);
/* wrbgf.c */
extern int write_bgf(FILE *file1, ums_type *mol);
/* wrbmin.c */
extern int write_bmin_com(FILE *file1, ums_type *mol);
/* wrc3d.c */
extern int write_chem3d2(FILE *file1, ums_type *mol);
extern int write_chem3d1(FILE *file1, ums_type *mol);
extern int write_mmads(FILE *file1, ums_type *mol);
extern int write_chem3d(FILE *file1, ums_type *mol, char *mol_typ);
/* wrcacao.c */
extern int write_caccrt(FILE *file1, ums_type *mol);
/* wrcache.c */
extern int write_cache(FILE *file1, ums_type *mol);
extern int SymbToNum(char *atSymb);
/* wrcacint.c */
extern int write_cacao_internal(FILE *file1, ums_type *mol);
extern void add_dummy_atoms(ums_type *mol);
extern void center_at_first_atom(ums_type *mol);
extern void set_hilderbrandt_connections(ums_type *mol);
extern void set_hilderbrandt_geometry(ums_type *mol);
/* wrchdrw.c */
extern int write_chem_draw(FILE *file1, ums_type *mol);
/* wrcontmp.c */
extern int write_conjure_tmplt(FILE *file1, ums_type *mol);
/* wrcsr.c */
extern int write_csr(FILE *fp, ums_type *mol);
extern void write_csr_header(FILE *fp, ums_type *mol);
extern void write_csr_coords(FILE *fp, ums_type *mol);
extern void write_size(int the_size, FILE *fp);
extern char *pad_string(char *input, int size);
/* wrcssr.c */
extern int write_cssr(FILE *file1, ums_type *mol);
extern coord_type calc_cell_dimensions(ums_type *mol);
extern void make_fractional(ums_type *mol, coord_type dim);
/* wrdock.c */
extern int write_dock(FILE *file1, ums_type *mol);
extern int total_heavy_atoms(ums_type *mol);
/* wrdpdb.c */
extern int write_dock_pdb(FILE *file1, ums_type *mol);
/* wrfeat.c */
extern int write_feat(FILE *file1, ums_type *mol);
/* wrfh.c */
extern int write_fenske_zmat(FILE *file1, ums_type *mol);
/* wrgamess.c */
extern int write_gamess_input(FILE *file1, ums_type *mol);
extern int write_gamess_cart(FILE *file1, ums_type *mol);
extern int write_gamess_zmatrix(FILE *file1, ums_type *mol);
extern int write_gamess_mopac(FILE *file1, ums_type *mol);
/* wrgau.c */
extern int write_gaussian(FILE *file1, ums_type *mol);
extern double fix_g90_angle(double angle);
extern int write_gaussian_template(FILE *file1, ums_type *mol);
extern void tmplt_cartgeom(ums_type *mol);
/* wrgaucrt.c */
extern int write_gaus_crt(FILE *file1, ums_type *mol);
/* wrhin.c */
extern int write_hin(FILE *file1, ums_type *mol);
/* wricon.c */
extern int write_icon8(FILE *file1, ums_type *mol);
/* wrint.c */
extern int write_mopac_internal(FILE *file1, ums_type *mol);
extern void cartint(ums_type *mol);
extern void cartgeom(ums_type *mol);
/* wrisis.c */
extern int write_isis(FILE *file1, ums_type *mol);
extern int structure_is_2D(ums_type *mol);
/* wrm3d.c */
extern int write_m3d(FILE *file1, ums_type *mol);
/* wrmaccs.c */
extern int write_maccs(FILE *file1, ums_type *mol);
/* wrmacmod.c */
extern int write_macromodel(FILE *file1, ums_type *mol);
extern int get_bond_order(ums_type *mol, int start, int end);
/* wrmcmol.c */
extern int write_mcmol(FILE *file1, ums_type *mol);
extern void get_atom_info(char *type_name, double *vdw_radius, double *bs_radius, double *red, double *grn, double *blu);
extern void translate_color(int *color, double *red, double *grn, double *blu);
/* wrmdl.c */
extern int write_molfile(FILE *file1, ums_type *mol);
/* wrmicro.c */
extern int write_micro(FILE *file1, ums_type *mol);
/* wrmimic.c */
extern int write_mac_mimic(FILE *file1, ums_type *mol);
/* wrmm2.c */
extern int write_mm2(FILE *file1, ums_type *mol);
extern int write_mouse(FILE *file1, ums_type *mol);
extern int write_mm2_input(FILE *file1, ums_type *mol);
extern int update_mm2_types(ums_type *mol, int i, int temp_type);
extern int check_for_carbonyl(ums_type *mol, int atom_num);
extern int check_for_cyclopropane(ums_type *mol, int atom1);
extern int bonded(ums_type *mol, int atom1, int atom2);
extern int type_mm2_hydrogens(ums_type *mol, int i);
/* wrmm3.c */
extern int write_mm3(FILE *file1, ums_type *mol);
/* wrmopac.c */
extern int write_mopac(FILE *file1, ums_type *mol);
/* wrpcmod.c */
extern int write_pcmod(FILE *file1, ums_type *mol);
/* wrpdb.c */
extern int write_pdb(FILE *file1, ums_type *mol);
extern int write_idatm(FILE *file1, ums_type *mol);
extern void assign_pdb_number(pdb_type_rec *pdb_types, int count);
/* wrpsgv.c */
extern int write_psgvb_cart(FILE *file1, ums_type *mol);
/* wrpsgvz.c */
extern int write_psgvb_zmat(FILE *file1, ums_type *mol);
/* wrsmiles.c */
extern int write_smiles(FILE *fp, ums_type *mol);
extern int ums_to_smiles_ct(ums_type *mol, smilescontab_t **ct);
extern void printcontab(smilescontab_t *contab);
/* wrspart.c */
extern int write_spartan(FILE *file1, ums_type *mol);
/* wrsybmol.c */
extern int write_sybyl_mol(FILE *file1, ums_type *mol);
/* wrsybyl.c */
extern int write_sybyl(FILE *file1, ums_type *mol);
extern void fix_sybyl_types(ums_type *mol);
/* wrtorlst.c */
extern int print_torsion_list(FILE *file1, ums_type *mol);
/* wrunichm.c */
extern int write_unichem(FILE *file1, ums_type *mol);
/* wrwiz.c */
extern int write_wizard(FILE *file1, ums_type *mol);
extern int compare_int(int *a, int *b);
/* wrxed.c */
extern int write_xed(FILE *file1, ums_type *mol);
/* wrxyz.c */
extern int write_xyz(FILE *file1, ums_type *mol);
/* wrmiv.c */
int write_mol_inventor(FILE *file1, ums_type *mol);

/* rdg96.c */
int gline(char *the_line, int nsz, FILE *f);
int cnt_atms(FILE *f);
int rd_sco_gr96(FILE *file1, ums_type *mol, double fac);
int read_gr96A(FILE *file1, ums_type *mol);
int read_gr96N(FILE *file1, ums_type *mol);
/* wrg96.c */
int wr_sco_gr96(FILE *file1, ums_type *mol, double fac);
int write_gr96A(FILE *file1, ums_type *mol);
int write_gr96N(FILE *file1, ums_type *mol);

/*wrtinker.c*/
int write_tinker(FILE *file1, ums_type *mol);

/* wrbox.c */
int write_box(FILE *fp, ums_type *mol);
void find_extents(ums_type *mol, vect_type *min, vect_type *max);

#ifdef __cplusplus
}
#endif

#endif  /* !__BABEL_BABEL_H__ */
