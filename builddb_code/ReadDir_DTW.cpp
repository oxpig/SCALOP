#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
//~ #include <thread>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <iomanip>
#include <sstream>

using namespace std;
using namespace arma;
using namespace boost::filesystem;


class Pdb_data {
    
    int anchor_len;
    int loop_len_1;
    int loop_len_2;
    int n_atoms;
    
    mat anchor_points_1;
    mat aligned_anchor_points_1;
    mat anchor_points_2;
    
    mat loop_points_1;
    mat aligned_loop_points_1;
    mat loop_points_2;
    
    rowvec anchor_atom_found_1;
    rowvec anchor_atom_found_2;
    
    rowvec residue_atom_found_1;
    rowvec residue_atom_found_2;
    
    void get_atom(char[], char[]);
    bool correct_atom(char[]);
    double get_double(string, int, int);
    int atom_index(char[]);
    int get_index(char*, char[][5], int);
    void read_number(char[], int, int, double*);
    void read_file(string, mat&, mat&, rowvec&, rowvec&, char[][5], int, int);
    void Kabsch(mat, mat, mat&, mat&, mat&);
    void Align_anchors();
    template <typename T> int sgn(T);
    void ResidueRM_TSD(mat, mat, rowvec, rowvec, double&, double&, int&);
    string format_number(double);
    
    
    public:
    
    Pdb_data(string, string, char[][5], int, int);
    void Retrieve_Residue(int, int, mat&, rowvec&);
    void calculate_DTW(double*);
    void write_modified_file(const char*, const char*, int);
    
};

void Pdb_data::get_atom(char line[], char atom[]){
    int atom_index = 0;
    for(int i=12; i<16; i++){
        if(isalnum(line[i])){
            atom[atom_index++] = line[i];
        }
    }
    atom[atom_index] = '\0';
}

bool Pdb_data::correct_atom(char atom[]){
    return strcmp(atom,"C")==0 || strcmp(atom,"CA")==0 || strcmp(atom,"N")==0 || strcmp(atom,"O")==0;
}

double Pdb_data::get_double(string line, int start, int length){
    double x = atof(line.substr(start, length).c_str());
    return x;
}

int Pdb_data::atom_index(char atom[]){
    if(strcmp(atom,"N")==0){
        return 0;
    } else if(strcmp(atom,"CA")==0){     
        return 1;
    } else if(strcmp(atom,"C")==0){
        return 2;
    } else if(strcmp(atom,"O")==0){
        return 3;
    }

    cout << "Unrecognized atom " << atom << endl;
}

int Pdb_data::get_index(char* atom, char atoms[][5], int n_atoms){
    for(int i=0; i<n_atoms; i++){
        if(strcmp(atom, atoms[i]) == 0){
            return i;
        }
    }
    return -1;
}

void Pdb_data::read_number(char line[], int from, int length, double *result){
    char num[length+1];
    memcpy(num, &line[from], length);
    num[length] = '\0';
    *result = atof(num);
}

void Pdb_data::read_file(string file, mat &anchor_points, mat &loop_points, rowvec &anchor_atom_found, rowvec &residue_atom_found, char atoms[][5], int n_atoms, int anchor_len){
    boost::filesystem::ifstream input(file);
    char buf[81];
    char c_aa[6], l_aa[6], atom[5];  
    char beginning[7];  
    const char *ATOM = "ATOM  ";
    const char *HETATM = "HETATM";
    int aa_index = -1;
    rowvec g_atom_found;
    mat points;
    
    c_aa[5] = '\0';
    l_aa[5] = '\0';
    atom[4] = '\0';
    beginning[4] = '\0';
    
    mat aa_atoms = zeros<mat>(3, n_atoms);
    //mat empty = zeros<mat>(n_atoms,3);
    rowvec l_atom_found = zeros<rowvec>(n_atoms);
    //rowvec empty_vec = zeros<rowvec>(n_atoms);
    
    // cout << "Reading file " << file << endl;
    
    while(input.good()){

        input.getline(buf, 81);
                
        memcpy(beginning, &buf[0], 6);

        beginning[6] = '\0';
        
        if (strcmp(beginning, ATOM) == 0 || strcmp(beginning, HETATM) == 0)
        {
        
            // Change Amino Acid
            memcpy(l_aa, &buf[22], 5);
            
            //cout << c_aa << endl;
            
            if(strcmp(c_aa, l_aa) != 0){
                //cout << "Change from" << c_aa << "to" << l_aa << " " << strcmp(c_aa, l_aa) << endl;
                memcpy(c_aa, &l_aa[0], 5);            
                aa_index++;
                //cout << aa_index << '\n';
                if(aa_index != 0){
                    points = join_horiz(points, aa_atoms);
                    
                    //cout << points << '\n';
                    
                    g_atom_found = join_horiz(g_atom_found, l_atom_found);
                }
                
                l_atom_found = zeros<rowvec>(n_atoms);
                aa_atoms = zeros<mat>(3, n_atoms);
                
            }
            
            get_atom(buf, atom);
            int atom_index = get_index(atom, atoms, n_atoms);
            
            if(atom_index != -1){
                double x, y, z;
                read_number(buf,30, 8,&x);
                read_number(buf,38, 8,&y);
                read_number(buf,46, 8,&z);
            
                colvec point = zeros<colvec>(3);
                point(0) = x;
                point(1) = y;
                point(2) = z;
            
            
                aa_atoms.col(atom_index) = point;
                
                l_atom_found(atom_index) = 1.0;
            }
            
        }
    }

    points = join_horiz(points, aa_atoms);
    
    g_atom_found = join_horiz(g_atom_found, l_atom_found);
    
    int points_len = points.n_cols;
    
    anchor_points = join_horiz(points.cols(0, anchor_len*n_atoms - 1), points.cols(points_len - anchor_len*n_atoms, points_len - 1));
    
    loop_points = points.cols(anchor_len*n_atoms, points_len - anchor_len*n_atoms - 1);
    
    anchor_atom_found = join_horiz(g_atom_found.subvec(0, anchor_len*n_atoms - 1), g_atom_found.subvec(points_len - anchor_len*n_atoms, points_len - 1));
    
    residue_atom_found = g_atom_found.subvec(anchor_len*n_atoms, points_len - anchor_len*n_atoms - 1);
    
    //~ cout << loop_points << endl;
    
}

void Pdb_data::Kabsch(mat P, mat Q, mat &anchor_mass_centre_1, mat &anchor_mass_centre_2, mat &Rotation_mat)
{
    mat Middle_mat;
    
    mat P_translated;
    mat Q_translated;
    mat P_rotated_translated;
    mat P_finish;
    
    mat U;
    vec s;
    mat V;
    
    mat R;
    
    colvec P_centroid;
    colvec Q_centroid;
    rowvec v_ones;
    
    mat Covariance_matrix;
    
    int d;
    int n_points;
    
    // Calculate the number of datapoints
    n_points = P.n_cols;
    
    //Calculate the coordinates of the centres of mass
    P_centroid=sum(P,1)/n_points;
    
    Q_centroid=sum(Q,1)/n_points;
    
    v_ones = ones<rowvec>(n_points);
    
    //Translate the points so that the cetre of mass coincides with the origin
    P_translated = P - P_centroid*v_ones;
    
    Q_translated = Q - Q_centroid*v_ones;
    
    //Calculate the covariance matrix   
    Covariance_matrix=P_translated*Q_translated.t()/n_points;
    
    //Calculate svd
    svd(U,s,V,Covariance_matrix);
    
    //Ensure a right-handed coordinate system
    d = sgn(det( V*U.t() ));
    
    Middle_mat = eye<mat>(3,3);
    
    Middle_mat(2,2) =  d;
    
    //Calculate the rotation matrix
    R=V*Middle_mat*U.t();
    
    //~ cout << Q << endl;
    //~ cout << R*(P - P_centroid*v_ones) + Q_centroid*v_ones << endl;
    
    aligned_anchor_points_1 = R*(P - P_centroid*v_ones) + Q_centroid*v_ones;
    
    anchor_mass_centre_1 = P_centroid;
    anchor_mass_centre_2 = Q_centroid;
    Rotation_mat = R;
    
}

void Pdb_data::Align_anchors(){
    
    mat anchor_mass_centre_1;
    mat anchor_mass_centre_2;
    mat Rotation_mat;
    int len_1;
    uvec logic_vec_1;
    
    Kabsch(anchor_points_1, anchor_points_2, anchor_mass_centre_1, anchor_mass_centre_2, Rotation_mat);
    
    rowvec v_ones = ones<rowvec>(sum(residue_atom_found_1));
    
    logic_vec_1 = find(residue_atom_found_1 == 1);
    
    aligned_loop_points_1 = loop_points_1;
    
    aligned_loop_points_1.cols(logic_vec_1) = aligned_loop_points_1.cols(logic_vec_1) - anchor_mass_centre_1*v_ones;
    
    aligned_loop_points_1.cols(logic_vec_1) = Rotation_mat*aligned_loop_points_1.cols(logic_vec_1);
    
    aligned_loop_points_1.cols(logic_vec_1) = aligned_loop_points_1.cols(logic_vec_1) + anchor_mass_centre_2*v_ones;
    
}

template <typename T> int Pdb_data::sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void Pdb_data::Retrieve_Residue(int data_index, int res_index, mat &Residue_points, rowvec &individual_atoms_found){
    
    if(data_index == 1)
    {
        Residue_points = aligned_loop_points_1.cols(res_index*n_atoms, (res_index + 1)*n_atoms - 1);
        
        individual_atoms_found = residue_atom_found_1.subvec(res_index*n_atoms, (res_index + 1)*n_atoms - 1);
        
    }
    else if(data_index == 2)
    {
        Residue_points = loop_points_2.cols(res_index*n_atoms, (res_index + 1)*n_atoms - 1);
        
        individual_atoms_found = residue_atom_found_2.subvec(res_index*n_atoms, (res_index + 1)*n_atoms - 1);
        
    }
    else
    {
        cout << "There are only two datasets" << endl;
        
    }
    
}

void Pdb_data::ResidueRM_TSD(mat Res1, mat Res2, rowvec individual_atoms_found_1, rowvec individual_atoms_found_2, double &RMSD, double &TSD, int &N){
    
    uvec logic_vec_both;
    
    mat Diff;
    
    rowvec total_individual_atoms_found;
    
    mat Res_inboth_1;
    
    mat Res_inboth_2;
    
    total_individual_atoms_found = individual_atoms_found_1 % individual_atoms_found_2;
    
    N = sum(total_individual_atoms_found);
    
    logic_vec_both = find(total_individual_atoms_found == 1);
    
    Res_inboth_1 = Res1.cols(logic_vec_both);
    
    Res_inboth_2 = Res2.cols(logic_vec_both);
    
    Diff = Res_inboth_1 - Res_inboth_2;
    
    RMSD = sqrt(1.0/N * sum(sum(Diff % Diff)));
    
    TSD = sum(sum(Diff % Diff));
    
}

void Pdb_data::calculate_DTW(double *result)
{
    
    double RMSD_val;
    
    double SD_val;
    
    int N;
    
    mat Residue_points_1;
    rowvec individual_atoms_found_1;
    
    mat Residue_points_2;
    rowvec individual_atoms_found_2;
    
    mat cost_mat = zeros<mat>(loop_len_1, loop_len_2);
    
    mat SD_mat = zeros<mat>(loop_len_1, loop_len_2);
    
    mat N_mat = zeros<mat>(loop_len_1, loop_len_2);
    
    uword index;
    
    double min_out_of_three;
    
    rowvec three_values = zeros<rowvec>(3);
    rowvec SD_three_values = zeros<rowvec>(3);
    rowvec N_three_values = zeros<rowvec>(3);
    
    Retrieve_Residue(1, 0 ,Residue_points_1, individual_atoms_found_1);
    Retrieve_Residue(2, 0 ,Residue_points_2, individual_atoms_found_2);
    
    ResidueRM_TSD(Residue_points_1, Residue_points_2, individual_atoms_found_1, individual_atoms_found_2, RMSD_val, SD_val, N);
    
    cost_mat(0,0) = RMSD_val;
    SD_mat(0,0) = SD_val;
    N_mat(0,0) = N;
    
    for(int i = 1; i < loop_len_1; i++)
    {
        Retrieve_Residue(1, i ,Residue_points_1, individual_atoms_found_1);
        
        ResidueRM_TSD(Residue_points_1, Residue_points_2, individual_atoms_found_1, individual_atoms_found_2, RMSD_val, SD_val, N);
        
        cost_mat(i,0) = cost_mat(i-1,0) + RMSD_val;
        
        SD_mat(i,0) = SD_mat(i-1,0) + SD_val;
        
        N_mat(i,0) = N_mat(i-1,0) + N;
        
    }
    
    Retrieve_Residue(1, 0 ,Residue_points_1, individual_atoms_found_1);
    
    for(int j = 1; j < loop_len_2; j++)
    {
        Retrieve_Residue(2, j ,Residue_points_2, individual_atoms_found_2);
        
        ResidueRM_TSD(Residue_points_1, Residue_points_2, individual_atoms_found_1, individual_atoms_found_2, RMSD_val, SD_val, N);
        
        cost_mat(0,j) = cost_mat(0, j-1) + RMSD_val;
        
        SD_mat(0,j) = SD_mat(0, j-1) + SD_val;
        
        N_mat(0,j) = N_mat(0, j-1) + N;
        
    }
    
    for(int i = 1; i < loop_len_1; i++)
    {
        // cout << "Here" << endl;
        for(int j = 1; j< loop_len_2; j++)
        {
            Retrieve_Residue(1, i ,Residue_points_1, individual_atoms_found_1);
            Retrieve_Residue(2, j ,Residue_points_2, individual_atoms_found_2);
            
            ResidueRM_TSD(Residue_points_1, Residue_points_2, individual_atoms_found_1, individual_atoms_found_2, RMSD_val, SD_val, N);
            
            three_values(0) = cost_mat(i-1, j);
            three_values(1) = cost_mat(i, j-1);
            three_values(2) = cost_mat(i-1, j-1);
            
            SD_three_values(0) = SD_mat(i-1, j);
            SD_three_values(1) = SD_mat(i, j-1);
            SD_three_values(2) = SD_mat(i-1, j-1);
            
            N_three_values(0) = N_mat(i-1, j);
            N_three_values(1) = N_mat(i, j-1);
            N_three_values(2) = N_mat(i-1, j-1);
            
            min_out_of_three = three_values.min(index);
            
            //cout << i << ' ' << j <<' '  << SD_val << three_values << endl;
            
            cost_mat(i,j) = min_out_of_three + RMSD_val;
            
            SD_mat(i, j) = SD_three_values(index) + SD_val;
            
            N_mat(i, j) = N_three_values(index) + N;
            
        }
    }
    
    //~ cout << cost_mat << endl;
    //~ 
    //~ cout << SD_mat << endl;
    //~ 
    //~ cout << N_mat << endl;
    
    N = N_mat(loop_len_1 - 1, loop_len_2 - 1);
    
    //~ cout << N << endl;
    //~ 
    //~ cout << sqrt(1.0/N * SD_mat(loop_len_1 - 1, loop_len_2 - 1)) << endl;
    
    *result = sqrt(1.0/N * SD_mat(loop_len_1 - 1, loop_len_2 - 1));
    
}

string Pdb_data::format_number(double Number_to_format)
{
    int wspc_missing;
    
    string Result;
    
    string complete_string;
    
    ostringstream Convert;
    
    Convert << fixed << setprecision(3) << Number_to_format;
    
    Result = Convert.str();
    
    wspc_missing = 7 - Result.length();
    
    complete_string = string(wspc_missing, ' ') + Result;
    
    return complete_string;

}

void Pdb_data::write_modified_file(const char *original_file, const char *result_file, int data_index)
{
    mat P_mat;
    
    if(data_index == 1)
    {
        P_mat = join_horiz(aligned_anchor_points_1.cols(0, anchor_len*n_atoms - 1), aligned_loop_points_1);
        P_mat = join_horiz(P_mat, aligned_anchor_points_1.cols(anchor_len*n_atoms, 2 * anchor_len*n_atoms - 1));
        
    }
    else if(data_index == 2)
    {
        P_mat = join_horiz(anchor_points_2.cols(0, anchor_len*n_atoms - 1), loop_points_2);
        P_mat = join_horiz(P_mat, anchor_points_2.cols(anchor_len*n_atoms, 2 * anchor_len*n_atoms - 1));
        
    }
    else
    {
        cout << "There are only two datasets" << endl;
        
    }
    
    
    boost::filesystem::ifstream infile(original_file);
    boost::filesystem::ofstream pdb_modified;
    
    double x,y,z;
    
    int i;
    
    pdb_modified.open(result_file);
    
    string line;
    
    string x_string;
    
    string y_string;
    
    string z_string;
    
    i=0;
    
    while(getline(infile, line))
    {
        
        if (line.compare(13, 3, "CA ")==0 or line.compare(13, 3, "N  ")==0 or line.compare(13, 3, "C  ")==0 or line.compare(13, 3, "O  ")==0)
        
        {
            x = P_mat(0,i);
            
            y = P_mat(1,i);
            
            z = P_mat(2,i);
            
            x_string = format_number(x);
            
            y_string = format_number(y);
            
            z_string = format_number(z);
            
            line.replace(31,7,x_string);
            
            line.replace(39,7,y_string);
            
            line.replace(47,7,z_string);
            
            //~ cout << line << line.length() << '\n';
            
            pdb_modified << line << '\n';
            
            i++;
        }
    }
    
    pdb_modified.close();
}
    
Pdb_data::Pdb_data(string file1, string file2, char atoms[][5], int n_atoms, int anchor_len){
    
    this->anchor_len = anchor_len;
    this->n_atoms = n_atoms;
    
    uvec logic_vec;
    rowvec total_anchor_atoms_found;
    
    
    //~ thread thread1(read_file, file1, ref(pos1), ref(atom_found1), atoms, n_atoms);
    //~ thread thread2(read_file, file2, ref(pos2), ref(atom_found2), atoms, n_atoms);
    //~ 
    //~ thread1.join();
    //~ thread2.join();

    // cout << file1 << endl;

    // cout << file2 << endl;
    
    read_file( file1, anchor_points_1, loop_points_1, anchor_atom_found_1, residue_atom_found_1, atoms, n_atoms, anchor_len);
    read_file( file2, anchor_points_2, loop_points_2, anchor_atom_found_2, residue_atom_found_2, atoms, n_atoms, anchor_len); 
    
    total_anchor_atoms_found = anchor_atom_found_1 % anchor_atom_found_2;
    
    // cout << anchor_points_1 << endl;
    
    // cout << anchor_points_2 << endl; 
    
    logic_vec = find(total_anchor_atoms_found == 1);
    
    anchor_points_1 = anchor_points_1.cols(logic_vec);
    
    anchor_points_2 = anchor_points_2.cols(logic_vec);
    
    //~ cout << anchor_points_1 << endl;
    
    if(anchor_points_1.n_elem != anchor_points_2.n_elem)
    {
        cout << "Anchor Lengths are not the same " << file1 << " " << file2 << endl;
    }
    
    //~ cout << anchor_points_2 << endl; 
    //~ 
    //~ cout << loop_points_1 << endl;
    //~ 
    //~ cout << loop_points_2 << endl;
    
    if(loop_points_1.n_cols % 4 != 0 || loop_points_2.n_cols % 4 != 0)
    {
        cout << "Something's wrong with the number of residues " << file1 << " " << file2 << endl;
    }
    
    loop_len_1 = loop_points_1.n_cols / 4.0;
    
    loop_len_2 = loop_points_2.n_cols / 4.0;
    
    //~ cout << loop_len_1 << " " << loop_len_2 << endl;
    
    Align_anchors();
    
    //~ cout << aligned_loop_points_1 << endl;
    //~ 
    //~ cout << loop_points_2 << endl;
    
    
    
    
}

int main(int argc, char* argv[])
{

    boost::filesystem::ofstream distmat_file;
    boost::filesystem::ofstream file_list;

    mat P1_kabsch;

    int anchor_len = 5;

    double result;
    
    char atoms[4][5];

    strcpy(atoms[0],"N");
    strcpy(atoms[1],"CA");
    strcpy(atoms[2],"C");
    strcpy(atoms[3],"O");

    string filename;

    string extension = "pdb";

    string current_extension;

    int filelen;

    int start;

    vector<string> filenames;

    mat distmat;

    int n_files;

    int row, col;

    int res_anchor_len = 5;
    int start_atom_anchor_len, end_atom_anchor_len;

    const char* directory = argv[1];

    const char* output_file_list_name = argv[2];

    const char* output_distmat_name = argv[3];

    path p(directory);

    // Read old file list
    vector<string> filenames_old;
    int n_oldfiles = 0;
    bool found_old = false;
    if (boost::filesystem::exists(output_file_list_name))
    {
        found_old = true;
	cout << "Found output file list\n";
        boost::filesystem::ifstream file_list_old(output_file_list_name);
        while(getline(file_list_old,filename))
        {
            filenames.push_back(filename);
        }
        n_oldfiles = filenames.size();
        cout << "Number of files: "<<n_oldfiles<<endl;
    }
    // Read old distmat
    mat distmat_old;
    if (found_old & boost::filesystem::exists(output_distmat_name))
    {
        boost::filesystem::ifstream distmat_file_old(output_distmat_name);
        distmat_old.load(distmat_file_old);
    }

    for (auto i = directory_iterator(p); i != directory_iterator(); i++)
    {
        if (!is_directory(i->path())) //we eliminate directories
        {

            filename = i->path().string();

            filelen = filename.length();

            start = filelen - 3;

            current_extension = filename.substr(start,filelen);

            if(!current_extension.compare(extension) & std::find(filenames.begin(),filenames.end(),filename) == filenames.end() ){
                filenames.push_back(filename);
            }

            
        }
        else
            continue;
    }

    cout << "Reading files and calculating distmat..." << endl;

    file_list.open(output_file_list_name);

    n_files = filenames.size();

    distmat = zeros<mat>(n_files, n_files);

    row = 0;

    for(auto i = filenames.begin(); i != filenames.end(); ++i)
    {

        //cout << row << " " << *i << endl;

        col = row + 1;
        file_list << *i << endl;

        for(auto j = i + 1; j != filenames.end(); ++j)
        {            
            //cout << *i << " " << *j << " " << row << " " << col << " ";
            if (found_old & col<n_oldfiles)
            {
		//cout << "skipped with "<<distmat_old(row,col)<<endl;
                distmat(row,col) = distmat_old(row,col);
                distmat(col,row) = distmat_old(col,row);
                col++;
                continue;
            }
            //cout << "calculating"<<endl;
            Pdb_data x(*i, *j, atoms, 4, anchor_len);

            x.calculate_DTW(&result);

            distmat(row,col) = result;

            distmat(col,row) = result;

            col++;

        }

        row++;
    }

    distmat_file.open(output_distmat_name);

    distmat_file << distmat;

    distmat_file.close();

    // cout << distmat << endl;



    // cout << "Done !" << endl;

    // Kabsch(all_Ps[0], all_Ps[2], &rmsd);

    // cout << filenames[0] << " " << filenames[1] << endl;

    // cout << rmsd << endl;
}
