// Pichenot Simon, Delifer Paul
/*
 * info602, TP2 : générateurs pseudo aléatoires simples
    Delifer Paul L3-Info*/

#include "tp2.h"
#include <assert.h>

/*
 * calcule l'inverse de a modulo m
 * (renvoie 0 si a n'est pas premier avec m ou si le résultat est incorrect)
 */
int64_t invert_mod(int64_t a, int64_t m)
{
    int64_t g = 0;
    int64_t x = 0;
    int64_t y = 0;
    gcd_bezout(&g, &x, &y, a, m);
    return (g == 1) ? mod(x, m) : 0; 
}

/*
 * craque un générateur congruentiel en cherchant le module, puis les nombres
 * a et c à partir de nombre générés
 *
 * Le tableau 'X' contient les 'nb' premiers nombres générés.
 * La fonction renvoie 1 si le générateur trouvé génère bien la liste donnée
 * et 0 sinon.
 */
int LCG_crack_ac(int nb, const int64_t *X, int64_t m, int64_t *a, int64_t *c)
{
    //On trouve les lignes suivantes après résolution à la main.
    // c = (X2-X1*(X3-X2))/(X2-X1) (x2-X1 !=0);
    // a = (X3-X2)/X2-X1 (X2-X1!=0);
    for (int i = 0; i < nb - 2; i++)
    {
        int invert = invert_mod(X[i + 1] - X[i], m);
        if (invert != 0)
        {
            *c = mod(X[i + 1] - X[i] * (X[i + 2] - X[i + 1]) * invert, m);
            *a = mod(((X[i + 2] - X[i + 1]) * invert), m);
            return 1;
        }
    }
    return 0;
}


/*
*Cette fonction permets de rendre plus lisible la fonction LCG_crack_m
* elle renvoie juste la différence entre Xi+1 et Xi.
*/
int iMoinsSuivant(int64_t const *X, int i)
{
    return X[i + 1] - X[i];
}

int LCG_crack_m(int nb, const int64_t *X, int64_t *m)
{
    int64_t a, c;
    if (nb < 5)
    {
        printf("\t # pas assez de valeur pour trouver m, le test échoue !\n");
        return 0;
    }
    else
    {
        printf("\t # OK, on peut calculer m !\n");
        *m = iMoinsSuivant(X, 2) * iMoinsSuivant(X, 0) - iMoinsSuivant(X, 1) * iMoinsSuivant(X, 1);
        int64_t tmp, tmp2;
        for (int i = 0; i < nb - 3; i++)
        {
            gcd_bezout(m, &tmp, &tmp2, *m, iMoinsSuivant(X, i + 2) * iMoinsSuivant(X, i) - iMoinsSuivant(X, i + 1) * iMoinsSuivant(X, i + 1));
        }
        LCG_crack_ac(nb, X, *m, &a, &c);
        if (LCG_crack_check(nb, X, *m, a, c))
        {
            return 1;
        }
        else
        {
            printf("\t# pas de solution pour générer ces valeurs\n");
            return 0;
        }
    }
}

int LCG_crack_check(int nb, const int64_t *X, int64_t m, int64_t a, int64_t c)
{
    // On vérifie qu'on peut générer toutes les valeurs de X avec les valeurs trouvées.
    for (int i = 0; i < nb - 1; i++)
    {
        // Si une valeur ne peut pas être générée alors on return 0
        if (mod(a * X[i] + c, m) != X[i + 1])
            return 0;
    }
    return 1;
}

/*
 * diagonalise une matrice sans changer la derniere colonne.
 *Le but de l'algorithme est de regarder si à l'indice i de la ligne i il y'a un 1.
 *Si oui on fait le xor avec les autres colonnes afin de supprimer les autres 1.
 *Cela revient a diagonaliser la matrice.
 */

int gauss(word *M, int nb_lignes)
{
    for (int i = 0; i < nb_lignes; i++)
    {
        if (!BIT(nb_lignes - i, M[i]))
        {
            for (int j = i; j < nb_lignes; j++)
            {
                if (BIT(nb_lignes - i, M[j])) // si on trouve un 1 on le remonte dans la diagonal.
                {
                    M[i] ^=  M[j];
                    break;
                }
                else if (j == nb_lignes - 1) // Si on ne trouve pas de 1 dans la colonne.
                {
                    return 0;
                }
            }
        }
        for (int j = 0; j < nb_lignes; j++)
        {
            if ((j != i) && BIT(nb_lignes - i, M[j]))
            {
                M[j] ^= M[i];
            }
        }
    }
    return 1;
}
/*
 * craque un générateur linéaire "fibonacci" en cherchant les "taps" qui
 * permettent de regénérer la suite.
 * Le tableau 'X' contient les 'nb' premiers bits générés par le
 * générateur.
 * La fonction renvoie 1 si les taps permettent de regénérer la suite, et 0
 * sinon.
 */
int LFSR_crack(int nb, const int *X, word *taps){
    word* diagX = (word*)X;
    int bit;
    for(int taille = 1; taille<=nb/2; taille++)
    {
        /
        if (gauss(diagX,taille))
        {
            for(int j = taille; j > 0; j--)
            {
                bit = BIT(taille,diagX[j]);   
                *taps ^= bit<<j;
            }
            return 1;
        }
        
    }
    return 0;
}