#include <bits/stdc++.h>
#define N 1100
#define M 1100
#define inf 2147483647
using namespace std;
int n,m,ans_for_LCSK,index_for_LCSK_in_B[N],length_for_LCSK_in_B[N];
int A[N], B[M], k, d1[N][M],d2[N][M],L1[N][M],L2[N][M],f[N][M],d[N][M],b1[N][M],b2[N][M];
//d[N][M] denote the length of the longest match between the prefixes of A[1 : i] = a1a2 ···ai and B[1 : j] = b1b2 ··· b j .

int Delta(int i)
{
    if(i>=k)
        return 1;
    return 0;
}
int Lambda(int i) 
{
    if(i-1<0)
        return 0;
    return (i - 1) % (2 * k) + 1;
}

void Epsilon(int i0,int i1,int j0,int j1)
{
    
    for (int i = 0; i <= k;i++)
        for (int j = 0; j <= m+k;j++)
            L1[i][j] = 0, d1[i][j] = 0,b1[i][j]=0;
    
    for (int i = i0; i <= i1;i++)
        for (int j = j0; j <= j1;j++)
        {
            
            L1[Lambda(i - 1)][j] = L1[0][j];
            if(A[i]==B[j])
                d1[Lambda(i)][j] = 1 + d1[Lambda(i - 1)][j - 1];
            L1[0][j] = L1[Lambda(i - 1)][j];
            if(d1[Lambda(i)][j]>=k)
            {
                b1[Lambda(i)][j] = 1 + b1[Lambda(i - 1)][j - 1];
                if(j-k>=j0&&i-k>=i0)
                    b1[Lambda(i)][j] = max(b1[Lambda(i)][j], k + L1[Lambda(i - k)][j - k]);
            }

            L1[0][j] = max(max(L1[Lambda(i - 1)][j], L1[0][j - 1]), b1[Lambda(i)][j]);
        }
    for (int j = j0; j <= j1;j++)
        L1[Lambda(i1)][j] = L1[0][j];

}
void Eta(int i0,int i1,int j0,int j1)
{   
    
    for (int i = 0; i <= k;i++)
        for (int j = 0; j <= m+k;j++)
            L2[i][j] = 0, d2[i][j] = 0,b2[i][j]=0;
    
    for (int i = i1; i >= i0;i--)
        for (int j = j1; j >= j0;j--)
        {
            L2[Lambda(i + 1)][j] = L2[0][j];
            if(A[i]==B[j])
                d2[Lambda(i)][j] = 1 + d2[Lambda(i + 1)][j + 1];
            L2[0][j] = max(L2[Lambda(i + 1)][j], L2[0][j + 1]);
            if(d2[Lambda(i)][j]>=k)
            {
                b2[Lambda(i)][j] = 1 + b2[Lambda(i + 1)][j + 1];
                if(j+k<=j1+1&&i+k<=i1+1)
                    b2[Lambda(i)][j] = max(b2[Lambda(i)][j], k + L2[Lambda(i + k)][j + k]);
            }

            L2[0][j] = max(max(L2[Lambda(i + 1)][j], L2[0][j + 1]), b2[Lambda(i)][j]);
        }
    for (int j = j0; j <= j1;j++)
        L2[Lambda(i0)][j] = L2[0][j];
}
void DandC(int i0,int i1,int j0,int j1)
{
    if(i1-i0+1<2*k)
        return;
    int l = max(2 * k - 1, (i1 - i0 + 1) / 2 + k);
    Epsilon(i0, i0 + l - 1, j0, j1);
    Eta(max(i0, i0 + l - 2 * k) + 1, i1, j0, j1);
    int s1 = inf, s2 = 0, tmp = 0;
    int k1, k2, l1, l2;
    // split(k1,k2,l1,l2, s1, s2,i1=i0+l-1, j0, j1)
    for (int i =  max(i0,i0 + l - 1 - 2*k + 1); i <= i0 + l - 1; i++)
        for (int j = j0 - 1; j <= j1;j++)
        {
            int t = L1[Lambda(i)][j] + L2[Lambda(i + 1)][j + 1];
            if(t>=tmp)
                tmp = t, k1 = i, k2 = j;
        }
    l1 = L1[Lambda(k1)][k2], l2 = L2[Lambda(k1 + 1)][k2 + 1];
    if(k<=l1&&l1<2*k)
    {
        for (int j = j0; j <= j1;j++)
        {
           
            if(L1[Lambda(k1)][j]==l1)
                s1 = min(s1, j);
        }
            
        s1 = s1 - l1 + 1;
    }
     if(k<=l2&&l2<2*k)
    {
       
        for (int j = j0; j <= j1;j++)
        {
            if(L2[Lambda(k1+1)][j]==l2)
                s2 = max(s2, j);
        }
            
    }
    if(l1>=2*k)
        DandC(i0, k1, j0, k2);
    else if(s1>0&&s1!=inf)
        index_for_LCSK_in_B[++ans_for_LCSK] = s1,length_for_LCSK_in_B[ans_for_LCSK]=l1;
    if(l2>=2*k)
        DandC(k1 + 1, i1, k2 + 1, j1);
    else if(s2>0)
        index_for_LCSK_in_B[++ans_for_LCSK] = s2,length_for_LCSK_in_B[ans_for_LCSK]=l2;
}

void LCSK()
{
    if(k>n/2)
        return;
    DandC(1, n, 1, m);
    printf("Answer for LCSK problem is: %d\n", ans_for_LCSK);
    printf("Substring begin with:");
    for (int i = 1; i <= ans_for_LCSK; i++)
        printf("%d ", index_for_LCSK_in_B[i]);
    printf("\nLength is: ");
    for (int i = 1; i <= ans_for_LCSK; i++)
        printf("%d ", length_for_LCSK_in_B[i]);
}
int main()
{
    printf("Input the length of A,B and k\n");
    scanf("%d%d%d", &n, &m, &k);
    printf("Input sequence A\n");
    for (int i = 1; i <= n;i++)
        scanf("%d", &A[i]);
    printf("Input sequence B\n");
    for (int i = 1; i <= m;i++)
        scanf("%d", &B[i]);
    LCSK();
    getchar();
    getchar();
    return 0;
}
