###################################################
###MÉTODOS PARA RESOLVER SISTEMAS DE ECUACIONES####
##################POR GABRIEL ROA##################
###################################################
###################CI: 25919459####################
#########Asignatura: Programación Numérica#########
############Profesor: Abelardo Monsalve############
####################UCLA - DCYT####################
###################################################

###
### MÉTODO DE GAUSS-JORDAN:
### Como eliminación gaussiana, pero resulta en la matriz identidad nxn en lugar de una triangular superior/inferior
### No necesita sustitución para resulver.
###

GaussJordan<-function(A,b) { ### 
	n=dim(A)                           

	if  (n[1]!=n[2]){
		return("Error: La matriz debe ser cuadrada")
	}

	d=det(A)

	if (d==0){
		return("Error: El sistema no tiene una única solución")
	}

	n=n[1]

	nb=n+1
	Ab=cbind(A,b)

	print("El sistema dado, representado de forma matricial, es el siguiente:")
	print(Ab)
	print("Apliquemos el método de Gauss-Jordan para resolver.")

	for (i in 1:n)
	{
		#Ab=PivoteoParcialDelante(Ab,i)  ### Descomentar para activar pivoteos.
		#Ab=PivoteoEscaladoDelante(Ab,i)
		for (j in 1:n){
			if (i!=j){
				factor=Ab[j,i]/Ab[i,i]
					if (factor!=0){
				 		Ab[j,]=Ab[j,]-factor*Ab[i,]
				 		print(Ab)		
				 	}
			 	}
			}
	 }

	print("Listo. Ahora vamos a resolver los cocientes restantes para hallar la solución del sistema:")
	x = numeric(n)

	for (i in 1:n){
		if (Ab[i,i]!=0){
			print(c("x",i,"=",Ab[i,nb],"/",Ab[i,i]),quote=FALSE)
			x[i]=Ab[i,nb]/Ab[i,i]
		}
	}

	print("La solución del sistema es: ")
	return(x)
}

###
### MÉTODO DE CRAMER:
### Incluye subrutina para resolver determinante de matriz 3x3,
### y programación general para resolver sistemas de ecuaciones 3x3 con el método de Cramer.
###

Determinante3x3<-function(A){
	a = A[1,1]*A[2,2]*A[3,3] ### Calcula los coeficientes positivos
	b = A[2,1]*A[3,2]*A[1,3]
	c = A[1,2]*A[2,3]*A[3,1]
	p = a+b+c

	a = A[3,1]*A[2,2]*A[1,3] ### Calcula los coeficientes negativos
	b = A[2,1]*A[1,2]*A[3,3]
	c = A[3,2]*A[2,3]*A[1,1]
	n = a+b+c

	return (p-n) ### Devuelve la diferencia
}

Cramer3x3<-function(A,b){ ### Calcula la solución del sistema de ecuaciones utilizando Regla de Cramer.
	n = dim(A)
	
	print("El sistema dado, representado de forma matricial, es el siguiente:")
	print(cbind(A,b))

	if (n[1] != n[2]){ ###Verifica que la matriz sea nxn
		return ("ERROR: La matriz debe ser cuadrada.")
	} else if (n[1] != 3){ ###Verifica que la matriz sea 3x3
		return ("ERROR: La matriz debe ser 3x3.")
	}

	print("Primero estudiemos el determinante de la matriz dada.")
	D = Determinante3x3(A)
	print(c("El determinante es igual a",D))

	if (D == 0){ ### Teorema: si el determinante de la matriz es igual a 0, entonces no tiene una única solución.
		print ("Atención: Como el determinante es 0, el sistema tiene infinitas soluciones o, en su defecto, no tiene.")
		print ("Para verificar, vamos a analizar los resultados de d1, d2 y d3.")
	}

	D1 = cbind(b,A[,2],A[,3]) ### Crea nuevas matrices para hallar sus determinantes
	D2 = cbind(A[,1],b,A[,3]) ### Los determinantes de estas matrices corresponden a los Di
	D3 = cbind(A[,1],A[,2],b)

	d1 = Determinante3x3(D1)
	d2 = Determinante3x3(D2)
	d3 = Determinante3x3(D3)

	if (D == 0 && d1 == 0 && d2 == 0 && d3 == 0){ ### En caso de que todos sean iguales a cero, entonces el sistema tiene infinitas soluciones
		print(c("d1 = ",d1))
		print(c("d2 = ",d2))
		print(c("d3 = ",d3))
		return("Como d1=d2=d3=0, podemos afirmar que el sistema tiene infinitas soluciones.")
	} else if (D == 0){ ### En caso de que todos sean distintos de cero, el sistema no tiene solución alguna.
		print(c("d1 = ",d1))
		print(c("d2 = ",d2))
		print(c("d3 = ",d3))
		return("Como d1 != d2 != d3, podemos afirmar que el sistema no tiene solución.")
	}

	x = c(d1/D,d2/D,d3/D) ### Caso contrario, si el sistema tiene soluciones, la calcula.

	print("Listo. La solución es:")
	print(x)
}

###
### MÉTODO DE CRAMER nxn:
### Método generalizado que resuelve la Regla de Cramer para sistemas de ecuaciones de nxn.
###

DeterminanteElimGauss<-function(A){ ### Cálculo del determinante de una matriz mediante aplicación de la eliminación gaussiana y cálculo de su traza.
	n=dim(A)[1]
	det=1

	for (i in n:2)
	{ ### Convertir matriz en triangular superior
		for (j in (i-1):1){
			factor=A[j,i]/A[i,i] 
			if (factor!=0){
		 		A[j,]=A[j,]-factor*A[i,]		
		 	}
		}
	}

	for (i in 1:n){ ### Producto de los elementos de la diagonal principal (traza)
	 	det = det*A[i,i]
	}

	return(round(det,5))
}

VerificarCerosCramer<-function(d){ ### Breve función que verifica que todos los elementos en un vector sean iguales a cero.
	for (i in 1 : length(d)){
		if (d[i]!=0)
			return(FALSE)
	}
	return(TRUE)
}

CramerNxN<-function(A,b){ ### Método de Cramer generalizado para matrices nxn
	n = dim(A)

	print("El sistema dado, representado de forma matricial, es el siguiente:")
	print(cbind(A,b))

	if (dim(A)[1] != dim(A)[2]){ ### Verifica que la matriz sea cuadrada
		return ("ERROR: La matriz debe ser cuadrada")
	}

	n = dim(A)[1]

	D = DeterminanteElimGauss(A) ### Determinante de la matriz dada
	print("Determinante de la matriz A del sistema: ")
	print(D)

	d = numeric(n)
	x = numeric(n)

	if (D == 0){ ### Si el determinante es igual a cero...
		print ("Atención: Como el determinante es 0, el sistema tiene infinitas soluciones o, en su defecto, no tiene.")
		print ("Para verificar, vamos a analizar los resultados de los distintos D.")
	}

	print("Comienza reemplazo para hallar los resultados de Di.")

	for (i in 1:n){ ### Reemplaza para las matrices D1, D2, ... , Dn
		aux = A
		aux[,i] = b ### Reemplaza la columna dada por la matriz columna b, que guarda los resultados
		print("Ahora trabajando con")
		print(c("D",i),quote=FALSE)
		print(aux)
		d[i] = round(DeterminanteElimGauss(aux),3) ### Halla el determinante y lo guarda en un arreglo d
		print("El determinante es igual a:")
		print(d[i])
	}

	if (VerificarCerosCramer(d)==TRUE && D == 0){ ### Si todo es igual a cero...
		print("Valores de D:")
		print(d)
		return("Como los valores de d son iguales a 0, podemos afirmar que el sistema tiene infinitas soluciones.")
	} else if (D == 0){ ### Si los valores de d son distintos de cero....
		print("Valores de D:")
		print(d)
		return("Como los valores de d son distintos, podemos afirmar que el sistema no tiene solución.")
	}

	print("Listos los reemplazos. Los valores de D son:")
	print(d)
	print("Ahora, vamos a efectuar la operación Di/D para hallar la solución, la cual es")

	for (i in 1:n) ### Calcula los distintos xi para las soluciones
		x[i] = d[i]/D

	print(round(x,3))
}

###
### ELIMINACIÓN GAUSSIANA: 
### Comenzamos con las rutinas de los pivoteos, y luego, el método de la eliminación gaussiana en cuestión.
###

### 
### ELIMINACIÓN GAUSSIANA HACIA ATRÁS
###

PivoteoParcialAtras<-function(Ab,i){ ### Funciona para eliminación hacia atrás
	n = dim(Ab)[1]
	g = abs(Ab[i,i])
	ind=i
	ind = which.max(abs(Ab[i,1:n])) ### Busca el índice del mayor módulo
	
	if (ind!=i){ ### Si el índice del mayor es distinto del índice del pivote, se hace cambio
		Faux=Ab[i,]
		Ab[i,]=Ab[ind,]
		Ab[ind,]=Faux[1:(n+1)]
		print("Intercambiando filas por Pivoteo Parcial: ")
		print(Ab)
	}

	return(Ab)
}

PivoteoEscaladoAtras<-function(Ab,i){ ### Pivoteo escalado para eliminación hacia atrás
	n=dim(Ab)[1]
	s=numeric(n)
	z=numeric(n)

	for (j in i : 1){ ### Guarda en vectores tanto los S de cada fila (mayor coeficiente) como los Z (valor del cociente)
		s[j]=max(abs(Ab[j,]))
		z[j]=abs(Ab[j,i]/s[j])
	}

	ind=which.max(z) ### Índice del mayor cociente escalado

	if (ind!=i){ ### Reemplazo si el índice de mayor no es el índice del pivote
		Faux=Ab[i,]
		Ab[i,]=Ab[ind,]
		Ab[ind,]=Faux[1:(n+1)]
		print("Intercambiando filas por Pivoteo Escalado: ")
		print(Ab)
	}

	return(Ab)
}

ElimGaussAtrasDelante<-function(A,b) { ### Eliminación Gaussiana: Eliminación hacia atrás, sustitución hacia delante
	n=dim(A)

	if  (n[1]!=n[2]){ ### Verifica que la matriz sea cuadrada
		return("Error: La matriz debe ser cuadrada")
	}

	d=det(A)

	if (d==0){ ### Verifica que el determinante sea distinto de 0
		return("Error: El sistema no tiene una única solución")
	}

	n=n[1] ### Guarda el tamaño de la matriz sin necesidad de acceder con un subíndice

	nb=n+1 ### Tamaño de la matriz ampliada
	Ab=cbind(A,b) ### Construye la matriz Ab

	print("El sistema dado, representado de forma matricial, es el siguiente:")
	print(Ab)
	print("Primero, vamos a convertir la matriz en una diagonal inferior...")

	for (i in n:2)
	{
		#Ab=PivoteoParcialAtras(Ab,i) ### Descomentar para activar el pivoteo
		#Ab=PivoteoEscaladoAtras(Ab,i)
		for (j in (i-1):1){
				factor=Ab[j,i]/Ab[i,i]
				if (factor!=0){
			 		Ab[j,]=Ab[j,]-factor*Ab[i,]
				 	print(Ab)								
			 	}
			}
	 }

	print("Listo. Ahora, vamos a sustituir para hallar la solución.")

	x=numeric(n)
	x[1]=Ab[1,nb]/Ab[1,1]

	for (i in 2 : n){
		suma=sum((Ab[i,1:(i-1)])*x[1:(i-1)])
		x[i]=(Ab[i,nb]-suma)/Ab[i,i]
	}
	
	print("¡Finalizado! La solución es:")
	print(x)
}


### 
### MÉTODO DE PIVOTEO ESCALONADO TRUNCADO A 3 DÍGITOS
###

trunc.dig <- function(x, digits) trunc(x*10^digits)/10^digits ### Crea método de truncamiento de n dígitos

PivoteoEscaladoAtrasTrunc3Dig<-function(Ab,i){ ### Pivoteo escalado para eliminación hacia atrás
	n=dim(Ab)[1]
	s=numeric(n)
	z=numeric(n)

	for (j in i : 1){ ### Guarda en vectores tanto los S de cada fila (mayor coeficiente) como los Z (valor del cociente)
		s[j]=trunc.dig(max(abs(Ab[j,])),3)
		z[j]=trunc.dig(abs(Ab[j,i]/s[j]),3)
	}

	ind=which.max(z) ### Índice del mayor cociente escalado

	if (ind!=i){ ### Reemplazo si el índice de mayor no es el índice del pivote
		Faux=Ab[i,]
		Ab[i,]=Ab[ind,]
		Ab[ind,]=Faux[1:(n+1)]
		print("Intercambiando filas por Pivoteo Escalado: ")
		print("Reemplaza la fila")
		print(Ab[i,])
		print("Por la fila")
		print(Ab[ind,])
		print("Resultado:")
		print(Ab)
	}

	return(Ab)
}

ElimGaussAtrasDelanteTrunc3Dig<-function(A,b) { ### Eliminación Gaussiana: Eliminación hacia atrás, sustitución hacia delante
	n=dim(A)

	if  (n[1]!=n[2]){ ### Verifica que la matriz sea cuadrada
		return("Error: La matriz debe ser cuadrada")
	}

	d=det(A)

	if (d==0){ ### Verifica que el determinante sea distinto de 0
		return("Error: El sistema no tiene una única solución")
	}

	n=n[1] ### Guarda el tamaño de la matriz sin necesidad de acceder con un subíndice

	nb=n+1 ### Tamaño de la matriz ampliada
	Ab=cbind(A,b) ### Construye la matriz Ab

	print("El sistema dado, representado de forma matricial, es el siguiente:")
	print(Ab)
	print("Primero, vamos a convertir la matriz en una diagonal inferior...")

	for (i in n:2)
	{
		Ab=trunc.dig(PivoteoEscaladoAtrasTrunc3Dig(Ab,i),3)
		for (j in (i-1):1){
				factor=trunc.dig(Ab[j,i]/Ab[i,i],3)
				print("Factor para volver 0:")
				print(factor)
				if (factor!=0){
			 		Ab[j,]=trunc.dig(Ab[j,]-factor*Ab[i,],3)
			 		print("La nueva matriz es:")
				 	print(Ab)								
			 	}
			}
	 }

	print("Listo. Ahora, vamos a sustituir para hallar la solución.")

	x=numeric(n)
	x[1]=trunc.dig(Ab[1,nb]/Ab[1,1],3)

	for (i in 2 : n){
		suma=trunc.dig(sum((Ab[i,1:(i-1)])*x[1:(i-1)]),3)
		x[i]=trunc.dig((Ab[i,nb]-suma)/Ab[i,i],3)
	}
	
	print("¡Finalizado! La solución es:")
	print(x)
}

### 
### MÉTODO DE PIVOTEO ESCALONADO REDONDEADO A 3 DÍGITOS
###

PivoteoEscaladoAtrasRnd3Dig<-function(Ab,i){ ### Pivoteo escalado para eliminación hacia atrás
	n=dim(Ab)[1]
	s=numeric(n)
	z=numeric(n)

	for (j in i : 1){ ### Guarda en vectores tanto los S de cada fila (mayor coeficiente) como los Z (valor del cociente)
		s[j]=round(max(abs(Ab[j,])),3)
		z[j]=round(abs(Ab[j,i]/s[j]),3)
	}

	ind=which.max(z) ### Índice del mayor cociente escalado

	if (ind!=i){ ### Reemplazo si el índice de mayor no es el índice del pivote
		Faux=Ab[i,]
		Ab[i,]=Ab[ind,]
		Ab[ind,]=Faux[1:(n+1)]
		print("Intercambiando filas por Pivoteo Escalado: ")
		print("Reemplaza la fila")
		print(Ab[i,])
		print("Por la fila")
		print(Ab[ind,])
		print("Resultado:")
		print(Ab)
	}

	return(Ab)
}

ElimGaussAtrasDelanteRnd3Dig<-function(A,b) { ### Eliminación Gaussiana: Eliminación hacia atrás, sustitución hacia delante
	n=dim(A)

	if  (n[1]!=n[2]){ ### Verifica que la matriz sea cuadrada
		return("Error: La matriz debe ser cuadrada")
	}

	d=det(A)

	if (d==0){ ### Verifica que el determinante sea distinto de 0
		return("Error: El sistema no tiene una única solución")
	}

	n=n[1] ### Guarda el tamaño de la matriz sin necesidad de acceder con un subíndice

	nb=n+1 ### Tamaño de la matriz ampliada
	Ab=cbind(A,b) ### Construye la matriz Ab

	print("El sistema dado, representado de forma matricial, es el siguiente:")
	print(Ab)
	print("Primero, vamos a convertir la matriz en una diagonal inferior...")

	for (i in n:2)
	{
		Ab=round(PivoteoEscaladoAtrasTrunc3Dig(Ab,i),3)
		for (j in (i-1):1){
				factor=round(Ab[j,i]/Ab[i,i],3)
				print("Factor para volver 0:")
				print(factor)
				if (factor!=0){
			 		Ab[j,]=round(Ab[j,]-factor*Ab[i,],3)
			 		print("La nueva matriz es:")
				 	print(Ab)								
			 	}
			}
	 }

	print("Listo. Ahora, vamos a sustituir para hallar la solución.")

	x=numeric(n)
	x[1]=round(Ab[1,nb]/Ab[1,1],3)

	for (i in 2 : n){
		suma=round(sum((Ab[i,1:(i-1)])*x[1:(i-1)]),3)
		x[i]=round((Ab[i,nb]-suma)/Ab[i,i],3)
	}
	
	print("¡Finalizado! La solución es:")
	print(x)
}

### 
### ELIMINACIÓN GAUSSIANA HACIA DELANTE
###

PivoteoParcialDelante<-function(Ab,i){ ### Pivoteo parcial para eliminación hacia delante.
	n=dim(Ab)[1]                       ### Funciona análogo al pivoteo parcial anterior, con cambios en los subíndices.
	g = abs(Ab[i,i])
	ind = which.max(abs(Ab[1:n,i]))

	if (ind!=i){
		Faux=Ab[i,]
		Ab[i,]=Ab[ind,]
		Ab[ind,]=Faux[1:(n+1)]
		print("Intercambiando filas por Pivoteo Parcial: ")
		print(Ab)
	}
	return(Ab)
}

PivoteoEscaladoDelante<-function(Ab,i){ ### Pivoteo escalado para eliminación hacia delante.
	n=dim(Ab)[1]                        ### Ídem.
	s=numeric(n)
	z=numeric(n)

	for (j in i : n){
		s[j]=max(abs(Ab[j,]))
		z[j]=abs(Ab[j,i]/s[j])
	}

	ind=which.max(z)

	if (ind!=i){
		Faux=Ab[i,]
		Ab[i,]=Ab[ind,]
		Ab[ind,]=Faux[1:(n+1)]
		print("Intercambiando filas por Pivoteo Escalado: ")
		print(Ab)
	}

	return(Ab)
}

ElimGaussDelanteAtras<-function(A,b) { ### Eliminación Gaussiana: Eliminación hacia delante, sustitución hacia atrás.
	n=dim(A)                           ### Nuevamente, análogo al caso anterior pero con cambios en los índices.

	if  (n[1]!=n[2]){
		return("Error: La matriz debe ser cuadrada")
	}

	d=det(A)

	if (d==0){
		return("Error: El sistema no tiene una única solución")
	}

	n=n[1]

	nb=n+1
	Ab=cbind(A,b)

	print("El sistema dado, representado de forma matricial, es el siguiente:")
	print(Ab)
	print("Primero, vamos a convertir la matriz en una diagonal superior...")

	for (i in 1:(n-1))
	{
		#Ab=PivoteoParcialDelante(Ab,i)  ### Descomentar para activar pivoteos.
		#Ab=PivoteoEscaladoDelante(Ab,i)
		for (j in (i+1):n){
				factor=Ab[j,i]/Ab[i,i]
				if (factor!=0){
			 		Ab[j,]=Ab[j,]-factor*Ab[i,]
			 		print(Ab)		
			 	}
			 	
			}
	 }

	print("Listo. Ahora, vamos a sustituir para hallar la solución, que viene dada por:")

	x=numeric(n)
	x[n]=Ab[n,nb]/Ab[n,n]

	for (i in (n-1) : 1){
		suma = sum((Ab[i,1:n])*x[1:n])
		x[i]=(Ab[i,nb]-suma)/Ab[i,i]
	}
	print(x)
}

###
### FACTORIZACIÓN LU:
### Método que realiza la factorización LU de una matriz y los posteriores cálculos para resolver el sistema.
###

ResolverLU<-function(A,b){ ### Método de resolución LU
	n=dim(A)

	if (n[1]!=n[2]){ ### Validaciones de siempre
		return("La matriz dada debe ser cuadrada.")
	}

	d=det(A)

	if (d==0){
		return("El sistema no tiene una única solución")
	}

	n=n[1]

	L=diag(1,n,n) ### Crea una matriz diagonal de orden nxn para L
	U=A           ### Copia la matriz original para la matriz U

	print("El sistema dado representado en forma matricial es:")
	print(cbind(A,b))
	cat("\n")
	print("===================================================")
	print("Primero realicemos la factorización A=L*U:")
	cat("\n")

	for (i in 1 : (n-1)){ ### Ciclo para factorización. Funciona igual a la eliminación hacia delante de la Elim. Gaussiana
		for (j in (i+1) : n){
			factor = U[j,i]/U[i,i]
			U[j,]=U[j,]-factor*U[i,]
			L[j,i]=factor ### Pero guarda el factor en la matriz L, para así construir la LU
			# print("L:")
			# print(L)
			# cat("\n")
			# print("U:")
			# print(U)
			# cat("\n")
		}
	}

	print("L: ")
	print(L)
	print("U: ")
	print(U)

	a=numeric(n)
	a[1]=b[1]

	print("===================================================")
	print("Ahora, resolvamos el subsistema asociado Ly=b")
	print(cbind(L,b))
	cat("\n")

	for (i in 2:n){ ### Sustitución hacia delante para resolver el subsistema de Ly=b
		sum=0
		for (j in 1:(i-1))
			sum=sum+L[i,j]*a[j]
		a[i]=b[i]-sum
	}

	print("La solución del subsistema es: ")
	print(a)
	cat("\n")
	print("===================================================")
	print("Ahora, resolvamos finalmente Ux=y")
	print(cbind(U,a))	
	cat("\n")

	x=numeric(n)
	x[n]=a[n]/U[n,n]
	
	for (i in (n-1) : 1){ ### Sustitución hacia atrás para resolver el sistema finalmente.
		sum = 0
		for (j in (i+1) : n){
			sum=sum+U[i,j]*x[j]
		}
		x[i]=(a[i]-sum)/U[i,i]
	}

	cat("\n")
	print("===================================================")
	print("La solución es: ")
	print("===================================================")
	print(x)
}

###
### MÉTODOS ITERATIVOS:
### Contiene los códigos con lo necesario para ejecutar Gauss-Seidel y Jacobi.
###

FuertementeDiagonal<-function(A,n){ ### Verifica que una matriz dada sea estrictamente dominante en sentido diagonal.
	for (i in 1:n){ ### Recorre toda la matriz
		sum = 0
		for (j in 1:n){
			if (i!=j){
				sum=sum+abs(A[i,j]) ### Recorre la matriz calculando la suma de los elementos de cada fila, ignorando la 
			}                       ### diagonal principal (i=j)
		}
		if (abs(A[i,i]<=sum)){ ### Verifica la condición
			print("La matriz no es estrictamente dominante en sentido diagonal.")
			return(0)
		}
	}
	print("La matriz es estrictamente dominante en sentido diagonal.")
	return(1) ### Retorna un valor booleano para poder imprimir un mensaje en el código más adelante
}

EvaluarError<-function(Xn,Xa,n){ ### Función que calcula el error para el criterio de parada en ambos métodos
	X=Xn-Xa ### Normaliza el vector para hallar la distancia
	dn=sqrt(sum(Xn[1:n]^2)) ### Halla la distancia nueva
	error=sqrt(sum(X[1:n]^2))/dn ### |Xn-Xa|/|Xn|
	return(error) ### Retorna el valor listo para realizar la comparación correspondiente.
}

###
### MÉTODO DE GAUSS-SEIDEL
### Iterativo, reemplaza los valores de X[i] a medida que los va obteniendo. Convergencia rápida.
###

GaussSeidel<-function(A,b,niter,tol){ ### Método de Gauss-Seidel. Recibe el número de iteraciones y la tolerancia.
	n = dim(A)

	if (n[1]!=n[2]){ ### Validaciones de todo sistema.
		return("La matriz debe ser cuadrada.")
	}

	n=n[1]
	d=det(A)

	if (d==0){
		return("El sistema no tiene una única solución")
	}

	print("El sistema dado es el siguiente:")
	print(cbind(A,b))
	cat("\n")

	if (FuertementeDiagonal(A,n)==0){ ### Verifica si el sistema _no_ es estrictamente diagonal.
		print("No se puede garantizar su convergencia.") ### Puede converger, pero lentamente, así como puede que no converja.
		cat("\n")
	}

	Xa=numeric(n) ### Crea un vector para almacenar los valores de x de la iteración anterior.
	Xn=numeric(n) ### Crea un vector para almacenar los valores de x de la iteración actual.

	for (i in 1 : n){
		Xn[i]=1
	}

	err=1 ### El primer error siempre es 1
	
	k=0 ### Índice para trabajar con la tolerancia.

	print("Comenzamos a iterar.")
	cat("\n")

	while((k<niter)&&(err>tol)){ ### Se hace mientras que no se alcance el máx de iteraciones y el error sea mayor que
		for (i in 1:n){			 ### la tolerancia.
			suma=0
			for (j in 1:n){
				if (i!=j){
					suma=suma+A[i,j]*Xn[j] ### Se podría usar un sum(), pero para ello habría que realizar la descomposición
				}                          ### Q-R. Notar que toma los valores del mismo vector de iteraciones nuevas.
			}                              ### Recordemos también que Xa y Xn fueron inicializados en 0.

			Xn[i]=(b[i]-suma)/A[i,i]       ### Calcula el nuevo X[i]
		}
		k=k+1 ### Incrementa el índice.
		print(cbind(("Iteración número:"),k)) 
		print("Valores de X:")
		print(Xn)
		cat("\n")
		err=EvaluarError(Xn,Xa,n) ### Realiza la evaluación del error.
		Xa=Xn                     ### Y después de eso, envía los resultados "nuevos" al vector de resultados anteriores.
	}

	if (k==niter){ ### Finalizó el ciclo, ¿por qué? Si k=niter, entonces se agotaron las iteraciones:
		print("Se ha alcanzado la cantidad máxima de iteraciones.")
		print("La solución hallada es la siguiente: ")
		print(Xn)
		cat("\n")
		print("Esta solución tiene un error relativo de: ")
		print(err)
	} else { ### Si no, entonces logré romper la tolerancia:
		print("Se ha conseguido una solución que satisface la tolerancia dada.")
		print("Esta solución es la siguiente: ")
		print(Xn)
		cat("\n")
		print("Su error relativo es: ")
		print(err)
	}
}

###
### MÉTODO DE JACOBI:
### Método iterativo que va evaluando de vector en vector en lugar de componente a componente.
###

Jacobi<-function(A,b,niter,tol){ ### Metodo de Jacobi. En su mayoría, análogo a Gauss-Seidel
	n = dim(A)

	if (n[1]!=n[2]){ ### Validaciones.
		return("La matriz debe ser cuadrada.")
	}

	n=n[1]
	d=det(A)

	if (d==0){
		return("El sistema no tiene una única solución")
	}

	print("El sistema dado es el siguiente:")
	print(cbind(A,b))
	cat("\n")

	if (FuertementeDiagonal(A,n)==0){ ### Misma idea que con Gauss-Seidel.
		print("No se puede garantizar su convergencia.")
		cat("\n")
	}

	Xa=numeric(n)
	Xn=numeric(n)

	for (i in 1 : n){
		Xa[i]=1
	}

	err=1
	
	k=0

	print("Comenzamos a iterar.")

	while((k<niter)&&(err>tol)){
		for (i in 1:n){
			suma=0
			for (j in 1:n){
				if (i!=j){
					suma=suma+A[i,j]*Xa[j] ### Acá el cambio: toma los valores de Xa (anteriores) y no de Xn (nuevos)
				}
			}
			Xn[i]=(b[i]-suma)/A[i,i]
		}
		k=k+1
		print(cbind(("Iteración número:"),k))
		print("Valores de X:")
		print(Xn)
		cat("\n")
		err=EvaluarError(Xn,Xa,n) ### Por lo demás, exactamente igual.
		Xa=Xn ### Sigue enviando Xn a Xa luego de evaluar el error.
	}

	if (k==niter){
		print("Se ha alcanzado la cantidad máxima de iteraciones.")
		print("La solución hallada es la siguiente: ")
		print(Xn)
		cat("\n")
		print("Esta solución tiene un error relativo de: ")
		print(err)
	} else {
		print("Se ha conseguido una solución que satisface la tolerancia dada.")
		print("Esta solución es la siguiente: ")
		print(Xn)
		cat("\n")
		print("Su error relativo es: ")
		print(err)
	}
}

###
### MÉTODO SOR
### Iterativo, reemplaza los valores de X[i] a medida que los va obteniendo. Trabajamos con un w para acelerar o garantizar la convergencia.
###

GaussSeidelSOR<-function(A,b,w,niter,tol){ ### Método de Gauss-Seidel. Recibe el número de iteraciones y la tolerancia.
	n = dim(A)

	if (n[1]!=n[2]){ ### Validaciones de todo sistema.
		return("La matriz debe ser cuadrada.")
	}

	n=n[1]
	d=det(A)

	if (d==0){
		return("El sistema no tiene una única solución")
	}

	print("El sistema dado es el siguiente:")
	print(cbind(A,b))
	cat("\n")

	if (FuertementeDiagonal(A,n)==0){ ### Verifica si el sistema _no_ es estrictamente diagonal.
		print("No se puede garantizar su convergencia.") ### Puede converger, pero lentamente, así como puede que no converja.
		cat("\n")
	}

	Xa=numeric(n) ### Crea un vector para almacenar los valores de x de la iteración anterior.
	Xn=numeric(n) ### Crea un vector para almacenar los valores de x de la iteración actual.

	for (i in 1 : n){
		Xn[i]=1
	}

	err=1 ### El primer error siempre es 1
	
	k=0 ### Índice para trabajar con la tolerancia.

	print("Comenzamos a iterar.")
	cat("\n")

	while((k<niter)&&(err>tol)){ ### Se hace mientras que no se alcance el máx de iteraciones y el error sea mayor que
		for (i in 1:n){			 ### la tolerancia.
			suma=0
			for (j in 1:n){
				if (i!=j){
					suma=suma+A[i,j]*Xn[j] ### Se podría usar un sum(), pero para ello habría que realizar la descomposición
				}                          ### Q-R. Notar que toma los valores del mismo vector de iteraciones nuevas.
			}                              ### Recordemos también que Xa y Xn fueron inicializados en 0.

			Xn[i]=(1-w)*Xa[i]+w*((b[i]-suma))/A[i,i]       ### Calcula el nuevo X[i] tomando en consideración la sobrerrelajación del w dado.
		}
		k=k+1 ### Incrementa el índice.
		print(cbind(("Iteración número:"),k)) 
		print("Valores de X:")
		print(Xn)
		cat("\n")
		err=EvaluarError(Xn,Xa,n) ### Realiza la evaluación del error.
		Xa=Xn                     ### Y después de eso, envía los resultados "nuevos" al vector de resultados anteriores.
	}

	if (k==niter){ ### Finalizó el ciclo, ¿por qué? Si k=niter, entonces se agotaron las iteraciones:
		print("Se ha alcanzado la cantidad máxima de iteraciones.")
		print("La solución hallada es la siguiente: ")
		print(Xn)
		cat("\n")
		print("Esta solución tiene un error relativo de: ")
		print(err)
	} else { ### Si no, entonces logré romper la tolerancia:
		print("Se ha conseguido una solución que satisface la tolerancia dada.")
		print("Esta solución es la siguiente: ")
		print(Xn)
		cat("\n")
		print("Su error relativo es: ")
		print(err)
	}
}
