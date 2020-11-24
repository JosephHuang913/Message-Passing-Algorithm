
void Encoder_H_Generator(int CodeType,int z_factor,int **H,int **inv_T);
void MatrixXORMultiply(int **A ,int **B,int **C,int rA,int cA,int cB);
void MS_decoder(int max_itr, int column, int row, double *receive_y, int **h, int *decoded_bit);
void SPA_decoder(int max_itr, int column, int row, double *receive_y, int **h, int *decoded_bit, double sigma_N);

