"""

    caixa(L, W)

Cria a forma de um envelope retangular de dimensões L x W.

    - 'L': comprimento do envelope retangular.
    - 'W': largura do envelope retangular.

Retorna um objeto do tipo Shape.

"""
caixa(L, W) = Shape([0, W, W, 0], [0, 0, L, L])


"""

    poligono(x_ref, y_ref, lista_poligono)

Cria a forma de um polígono referenciado a partir do ponto (x_ref, y_ref).

    - 'x_ref': abscissa do ponto de referência.
    - 'y_ref': ordenada do ponto de referência.
    - 'lista_poligono': matrix Px2, cuja primeira coluna representa a coordenada `x` e a segunda coluna a coordenada `y` de cada ponto do polígono

Retorna um objeto do tipo Shape.

"""
poligono(x_ref, y_ref, lista_poligono) = Shape(x_ref .+ lista_poligono[:, 1], 
    y_ref .+ lista_poligono[:, 2])

poligono(lista_poligono) = poligono(0.0, 0.0, lista_poligono)

function uniao_poligonos(P)

    return [poligono(P[i]) for i = 1:length(P)]

end

"""

    poligono_rotacionado(x_ref, y_ref, θ, lista_poligono)

Cria a forma de um polígono referenciado a partir do ponto (x_ref, y_ref) e 
rotacionado com ângulo θ.

    - 'x_ref': abscissa do ponto de referência.
    - 'y_ref': ordenada do ponto de referência.
    - 'θ': ângulo de rotação.
    - 'lista_poligono': lista de vértices do polígono.

Retorna um objeto do tipo Shape.

"""
poligono_rotacionado(x_ref, y_ref, θ, lista_poligono) = Shape(
    x_ref .+ cos(θ) * lista_poligono[:, 1] .- sin(θ) * lista_poligono[:, 2], 
    y_ref .+ sin(θ) * lista_poligono[:, 1] .+ cos(θ) * lista_poligono[:, 2])


quadrado(v, l) = poligono_rotacionado(v[1], v[2], 0.0, vcat([-l/2 -l/2], [l/2 -l/2], [l/2 l/2], [-l/2 l/2]))
quadrado(v, l, θ) = poligono_rotacionado(v[1], v[2], θ, vcat([-l/2 -l/2], [l/2 -l/2], [l/2 l/2], [-l/2 l/2]))

_triangulo_base = let 

    # baricentro do triangulo com vértice em (0, 0)
    bc = [1/2 sqrt(3) / 6]

    [
        0.0  0.0;
        1.0  0.0;
        1/2 (1 * sqrt(3) / 2)
    ] .- bc
    
end

triangulo_eq(d, θ, l) = triangulo(d, θ, l * _triangulo_base)

triangulo(d, θ, triangulo_base) = poligono_rotacionado(d[1], d[2], θ, triangulo_base)

"""

    circulo(x_ref, y_ref, raio)

Cria a forma de um círculo referenciado a partir do ponto (x_ref, y_ref) e de raio 
'raio'.

    - 'x_ref': abscissa do centro de referência.
    - 'y_ref': ordenada do centro de referência.
    - 'raio': raio do círculo.

Retorna um objeto do tipo Shape.

"""
function circulo(x_ref, y_ref, raio)
    θ = 0:5:360
    Shape(raio * sind.(θ) .+ x_ref, raio * cosd.(θ) .+ y_ref)
end

circulo(c, r) = circulo(c[1], c[2], r)

# ------------------------------------------------------------------------------
# Funções para desenhar problema e solução.
# ------------------------------------------------------------------------------


"""

    translacao_poligono(lista_vertices_pol)

Translada cada polígono até a origem.

    - 'lista_vertice_pol': lista contendo a lista de vértices de cada polígono.

"""
function translacao_poligono(lista_vertices_pol)

    # Faz uma cópia da lista de vértices
    nova_lista = deepcopy(lista_vertices_pol)
    
    # Calcula o número de polígonos
    np = length(nova_lista)
    
    for i = 1:np

        # Calcula o número de componentes convexas do polígono i
        ncc = length(nova_lista[i])

        for j = 1:ncc

            # Calcula o número de vértices da componente convexa j do polígono i
            nv = size(nova_lista[i][j])[1]

            # Subtrai as coordenadas do primeiro vértice da primeira componente 
            # convexa do polígono i dos demais vértices.
            for k=1:nv

                nova_lista[i][j][k, :] -= lista_vertices_pol[i][1][1, :]

            end

        end
    
    end
    
    return nova_lista

end


"""

    plota_problema(lista_centros_circ, lista_raios_circ, lista_vertice_pol, L, W)

Plota o problema de corte e empacotamento.

    - 'lista_centros_circ': lista contendo a lista dos centros dos círculos.
    - 'lista_raios_circ': lista contendo os raios dos círculos.
    - 'lista_vertice_pol': lista contendo a lista de vértices de cada polígono.
    - 'L': comprimento do envelope retangular.
    - 'W': largura do envelope retangular.

"""
function plota_problema(lista_vertice_pol, L, W)

   # nc = size(lista_raios_circ)[1]
    np = length(lista_vertice_pol)
    print(np)
    # Plota o envelope retangular.
    fig = plot(caixa(L, W), fillcolor=:white, legend = false, aspect_ratio=:equal)

    # Plota os círculos.
    #for i = 1:nc
       # plot!(fig, circulo(lista_centros_circ[i, 1], lista_centros_circ[i, 2], lista_raios_circ[i]), fillcolor = plot_color(:red, 0.5))
    #end

    # Plota os polígonos.
    for i = 1:np
        # Calcula o número de componentes convexas do polígono i.
        ncc = length(lista_vertice_pol[i])

        # Plota cada componente convexa.        
        for j = 1:ncc
            plot!(fig, poligono(0.0, 0.0, lista_vertice_pol[i][j]), fillcolor = plot_color(:blue, 0.5))
        end
    end

    fig

end


function reta(lista_poligonos)
    np = length(lista_poligonos)
    rt = Vector{Any}(undef,np)
    for i = 1:np
        ncc = length(lista_poligonos[i])
        rt[i]= Vector{Any}(undef,ncc)
        for j = 1 :ncc
        rt[i][j] = lista_poligonos[i][j][1:2, :]
        end
    end
         return rt
end

    
"""

    plota_solucao(lista_raios_circ, lista_vertice_pol, L, W, x)

Plota a solução do problema de corte e empacotamento.

    - 'lista_raios_circ': lista contendo os raios dos círculos.
    - 'lista_vertice_pol': lista contendo a lista de vértices de cada polígono (translado para origem).
    - 'L': comprimento do envelope retangular.
    - 'W': largura do envelope retangular.
    - 'x': solução encontrada pelo solver no formato [xC, yC, xP, yP, θP].

"""
function plota_solucao(lista_raios_circ, lista_vertice_pol, L, W, x)

    nc = size(lista_raios_circ)[1]
    np = length(lista_vertice_pol)

    # Plota o envelope retangular.
    fig = plot(caixa(L, W), fillcolor = :white, legend = false, aspect_ratio = :equal)

    # Plota os círculos.
    for i = 1:nc
        plot!(fig, circulo(x[2*i-1], x[2*i], lista_raios_circ[i]), fillcolor = plot_color(:red, 0.5))
    end

    # Plota os polígonos.
    for i = 1:np
        # Deslocamento no vetor de soluções 'x'.
        k = 2 * nc + 3 * i

        # Calcula o número de componentes convexas do polígono i.
        ncc = length(lista_vertice_pol[i])

        # Plota cada componente convexa.
        for j = 1:ncc
            plot!(fig, poligono_rotacionado(x[k - 2], x[k - 1], x[k], lista_vertice_pol[i][j]), fillcolor = plot_color(:blue, 0.5))
        end
    end

    fig

end

desenha_solucao_cq(raios, xC, yC, W) = desenha_solucao_cq(raios, xC, yC, W, W)

function desenha_solucao_cq(raios, xC, yC, W, L)

    @assert length(raios) == length(xC) == length(yC)

    fig = plot(caixa(L, W), c=:dodgerblue, lc=:black, lw=3, label=false, aspect_ratio=true)

    for (c, r) in zip(zip(xC, yC), raios)

        plot!(fig, circulo(c, r), c=false, lc=:firebrick, lw=3, label=false, aspect_ratio=true)

    end

    return fig

end

function desenha_solucao_qc(ℓ, xQ, yQ, θ, R)

    @assert length(xQ) == length(yQ) == length(θ)

    fig = plot(circulo([0, 0], R), c=:dodgerblue, lc=:black, lw=3, label=false, aspect_ratio=true)

    for (v, a) in zip(zip(xQ, yQ), θ)

        plot!(fig, quadrado(v, ℓ, a), c=false, lc=:firebrick, lw=3, label=false, aspect_ratio=true)

    end

    return fig

end