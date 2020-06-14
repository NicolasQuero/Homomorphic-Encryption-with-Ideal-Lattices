# Homomorphic-Encryption-with-Ideal-Lattices

Dans le cadre d'un projet pour le M1 MIC à Paris Diderot, nous avons cherché à implémenter la version fully homomorphic proposée par Craig Gentry dans son papier "Implementing Gentry's Fully-Homomorphic Encryption Scheme" publié en 2011.

Après un long travail de lecture et d'apprentissage sur le chiffrement homomorphe, nous avons pu créer le système de chiffrement SWHE par les réseaux idéaux qui permet alors d'effectuer un nombre limité d'opérations sur les chiffrés. Mais pour passer ce système en fully homomorphic, il nous reste encore à implémenter les dernières étapes du squashing permettant d'alléger la méthode de déchiffrement et de permettre de passer des indices sur la clé secrète pour effectuer le bootstrapping.
