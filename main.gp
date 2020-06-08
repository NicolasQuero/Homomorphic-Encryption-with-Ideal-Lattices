gen_ideal(n, t) = {
	my(f=0);
	for(i=1,n,
		f = f + (random(1 << t) - (1 << (t-1))) * X^(i-1);
		\\print(f);
	);
	f;
}

vec_to_pol(v, n) = {
	my(f = 0);
	for(i=1, n+1,
		f = f + v[i]*(X^(i-1));

	);
	f;
}

pol_to_vec(f, n) = {
	my(v = vector(n+1));
	for(i=1, n+1,
		v[i] = polcoef(f, i-1);
	);
	v;
}

secret_key(g, f, n) = { \\ n la dimension du reseau, g le gen de l'ideal, f le polynome de quotient
	my(m);
	m = matrix(n, n, i, j, polcoef(lift(Mod((X^(i-1))*g, f)), j-1));
	m;
}

random_vec(n) = {
	my(v = vector(n));
	for(i=1, n,
		v[i] = random(2);
	);
	v;
}

mod_lattice(v, B) = {
	my(e);
	e = v - B*round((B^-1)*v);
	e;
}

encrypt(b, pk, n) = {
	my(r);
	my(e);
	e = vector(n);
	r = random_vec(n);
	e = 2*r;
	e[1] = e[1] + b;
	mod_lattice(e~, pk);
}

decrypt(c, w, d) = {
	res = lift(Mod(c[1], d)*w);
	if(res >= d/2, res -= d);
	res%2;
}

mult_vec(v1, v2, n) = {
	return(concat(v1[1]*v2[1],v1[2..n])); \\ On considère que ce sont des chiffrés
}

main() = {
	n = 2^6;
	t = 2^6;
	S = 2^5;
	s = 15; \\ poids du vecteur sig
	d = 0; \\ On génère l'idéal, la base secrète et un d nécessairement impair
	while(d%2 == 0,
		g = gen_ideal(n, t);
		f = X^n + 1;
		sk = secret_key(g, f, n);
		d = matdet(sk);
	);

	v = sk[1,]; \\ v est la première ligne de sk qui se décale à chaque ligne
	vp = vec_to_pol(v, n-1);
	xgcd = gcdext(vp, f)*d;
	wx = pol_to_vec(xgcd[1], n-1); \\ wx est le polynôme inverse de v retrouvé dans W
	w = 1;
	for(i = 1, n,
		if(wx[i]%2 == 1 && wx[i] > 1, w = wx[i]; iw = i; break;);
	); \\ On trouve un w impair qui est alors la clé privée

	r = (wx[1]*(wx[2]^-1))%d; \\ paramètre dans pk
	pk = mathnf(sk); \\ clé publique
	R = random(d); while(gcd(R, d) != 1, R = random(d));

	sig = generate_sigma(s, S);
	sig2 = generate_sig2(s, S); \\ optimisé
	\\ génère les ensembles xk et sig tels que sum(x(i)*sig(i)) = w
	xk = generate_xk(w, d, sig, s, S);
	xk2 = generate_xk2(w, d, sig2, R, s, S); \\ optimisé

	/* Test de 10 décryptages aléatoires */

	for(i = 1, 10,
		a = random(2);
		c = encrypt(a, pk, n);
		m = decryption_squashed(c, d, xk, sig, S)%2;
		if(m != a, print("Problème de décryptage"), print("Décryptage de ", a, " avec succès"));

	);
}

generate_sigma(s, S) = {
	sig = vector(S);
	my(index);
	for (i = 1, s,
		index = random(S-1)+1;
		while(sig[index] != 0, index = random(S-1)+1);
		sig[index] = 1;
	);
	sig;
}

generate_sig2(s, S) = { \\ amélioration du squashing
	sig2 = matrix(s, S);
	my(used_indexes, row, index);
	used_indexes = [];
	row = 2;
	while(#used_indexes < s-1,
		index = random(S-2) + 2;
		if(vecsearch(used_indexes, index) == 0,
			sig2[row, index] = 1;
			used_indexes = concat(used_indexes, index);
			vecsort(used_indexes);
			row++;
		);
	);
	sig2[1,1] = 1;
	sig2;
}

generate_xk(w, d, sig, s, S) = {
	xk = vector(S, x, random(d));
	my(compteur, i);
	sum_xk = Mod(0, d);
	compteur = 0;
	i = 1;
	while(compteur < s-1,
		if(sig[i] == 1, sum_xk += xk[i]; compteur++;);
		i++;
	);
	while(i <= S,
		if(sig[i] == 1,
			xk[i] = lift(w - sum_xk);
			sum_xk += xk[i];
			compteur++;
		);
		i++;
	);
	return(xk);
}

generate_xk2(w, d, sig, R, s, S) = { \\ amélioration du squashing
	xk2 = vector(s);
	my(row, sum_xk2);
	row = 2;
	sum_xk2 = Mod(0, d);
	for(i = 2, s,
		j = vecsearch(sig, 1);
		xk2[i] = random(d);
		sum_xk2 += xk2[i]*(Mod(R, d)^j);
	); \\ xk2[1]*R + sum_xk2 = w
	xk2[1] = (w - sum_xk2)*Mod(R, d)^-1;
	sum_xk2 += xk2[1];
	xk2;

}

decryption_squashed(c, d, xk, sig, S) = {
	my(dec);
	dec = Mod(0, d);
	for (i = 1, S, dec += sig[i]*Mod(c[1]*xk[i], d);); \\ sum(x(i)*sig(i)) = w
	dec = lift(dec);
	while (dec > d>>2, dec -= d);
	dec;
}

binary_round(x, p) = { \\ 2 > x >= 0, p precision parameter
	res = [];
	while(#res < p+1,
		if (x >= 1,
			x-= 1;
			res = concat(res, 1);
			x *= 2;,
			\\ else
			res = concat(res, 0);
			x *= 2;
		);
	);
	res;
}

test_operations(nb) = {
	a = random(2);
	expected = a;
	c = encrypt(a, pk, n);
	for(i = 1, nb,
		a = random(2);
		b = random(2);
		if(b==1, c += encrypt(a, pk, n); expected = (expected + a)%2; print(" + ", a);,
		 c = mult_vec(c, encrypt(a, pk, n), n); expected = expected*a; print(" * ", a););
		m = decryption_squashed(c, d, xk, sig, S)%2;
		if(m == expected, print(i, " : Attendu ", expected, ", Obtenu : ", m),
		 print(i, " : Attendu ", a, ", Obtenu : ", m, " (Problème)"); break;);
	);
}

/*decryption_squashed2(c, d, xk, sig, R) = { \\ en construction
	my(dec, y, z, p);
	size = matsize(sig);
	s = size[1];
	S = size[2];
	y = matrix(s, S);
	z = matrix(s, S);
	p = ceil(log(s+1)/log(2));
	dec = Mod(0, 2);

	for(k = 1, s, \\ calcul des y_k,i
		for(i = 1, S,
			y[k, i] = lift(Mod(c[1]*xk[k], d)*R^i);
			z[k, i] = binary_round(y[k, i], p);
		);
	);

	my(mult);
	mult = matrix(s, S);
	for(i = 1, s,
		for(j=1, S,
			if(z[i,j] != 0,
				\\ effectuer la multiplication sigma * z
			);
		);
	);
}

gsa(var) = {
	size = matsize(var);
	s = size[1];
	p = size[2];
	for(row = 1, s,
		stack = 0;
		for(col = 1, p,

		);
	);
}
*/
