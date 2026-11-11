# SciCalc Pro — motor simbólico isolado (Pyodide/SymPy)
# © 2026 João Carlos Barumby. Integração adicional ao SciCalc.
# O navegador executa este arquivo num Web Worker; nenhum cálculo vai a um servidor.
import json
import re
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, convert_xor, implicit_multiplication_application
from sympy.calculus.util import continuous_domain

_X = sp.Symbol('x', real=True)
_FUNCS = {
    'x': _X, 'pi': sp.pi, 'E': sp.E, 'e': sp.E,
    'oo': sp.oo, 'inf': sp.oo,
    'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan,
    'asin': sp.asin, 'acos': sp.acos, 'atan': sp.atan,
    'sinh': sp.sinh, 'cosh': sp.cosh, 'tanh': sp.tanh,
    'exp': sp.exp, 'log': sp.log, 'ln': sp.log, 'sqrt': sp.sqrt,
    'abs': sp.Abs, 'Abs': sp.Abs, 'erf': sp.erf,
    'gamma': sp.gamma,
}
_GLOBALS = {
    '__builtins__': {}, 'Integer': sp.Integer, 'Float': sp.Float,
    'Rational': sp.Rational, 'Symbol': sp.Symbol,
    'Add': sp.Add, 'Mul': sp.Mul, 'Pow': sp.Pow,
}
_TRANSFORMS = standard_transformations + (convert_xor, implicit_multiplication_application)
_IDENTIFIERS = set(_FUNCS)

def _parse(value, label='expressão'):
    source = str(value).strip()
    if not source or len(source) > 240:
        raise ValueError(f'{label}: digite até 240 caracteres.')
    if not re.fullmatch(r'[A-Za-z0-9+\-*/^()., \t]+', source):
        raise ValueError(f'{label}: use apenas números, x e funções matemáticas permitidas.')
    tokens = re.findall(r'[A-Za-z][A-Za-z0-9]*', source)
    if any(token not in _IDENTIFIERS for token in tokens):
        unknown = sorted(set(token for token in tokens if token not in _IDENTIFIERS))
        raise ValueError('Identificador não permitido: ' + ', '.join(unknown))
    if re.search(r'(\*\*|\^)\s*[0-9]{4,}', source):
        raise ValueError('Expoente excessivamente grande.')
    try:
        result = parse_expr(source, local_dict=_FUNCS, global_dict=_GLOBALS,
                            transformations=_TRANSFORMS, evaluate=True)
    except Exception as exc:
        raise ValueError(f'Não foi possível interpretar {label}: {type(exc).__name__}.')
    if not isinstance(result, sp.Basic):
        raise ValueError('A entrada não representa uma expressão matemática.')
    if result.free_symbols - {_X}:
        raise ValueError('Esta versão aceita apenas a variável x.')
    return result

def _tex(value):
    return sp.latex(value)

def _numeric(value):
    if value in (sp.oo, -sp.oo, sp.zoo, sp.nan):
        return None
    try:
        estimate = sp.N(value, 16)
        if estimate.is_real and estimate.is_finite:
            return str(estimate)
    except Exception:
        pass
    return None

def _result(expr, answer, status, note, **extra):
    output = {
        'ok': True, 'input_tex': _tex(expr),
        'result_tex': _tex(answer), 'result_text': str(answer),
        'approx': _numeric(answer), 'status': status, 'note': note
    }
    output.update(extra)
    return output

def scicalc_pro_solve(request_json):
    req = json.loads(request_json)
    action = str(req.get('action', ''))
    expr = _parse(req.get('expr', ''), 'função')
    if action == 'derivative':
        order = int(req.get('order', 1))
        if not 1 <= order <= 6:
            raise ValueError('Escolha uma ordem de derivação entre 1 e 6.')
        value = sp.diff(expr, _X, order)
        output = _result(expr, value, 'Simbólico', f'Derivada de ordem {order}, calculada por diferenciação simbólica.')
        point = str(req.get('point', '')).strip()
        if point:
            p = _parse(point, 'ponto de avaliação')
            if p.free_symbols or p.is_real is not True:
                raise ValueError('O ponto deve ser um valor real.')
            evaluated = value.subs(_X, p)
            output['point_tex'] = _tex(p)
            output['point_result_tex'] = _tex(evaluated)
            output['point_approx'] = _numeric(evaluated)
        return json.dumps(output, ensure_ascii=False)

    if action == 'integral':
        low = str(req.get('lower', '')).strip()
        high = str(req.get('upper', '')).strip()
        if bool(low) != bool(high):
            raise ValueError('Informe ambos os limites ou deixe os dois em branco.')
        if not low:
            answer = sp.integrate(expr, _X)
            if answer.has(sp.Integral):
                output = _result(expr, answer, 'Não resolvida em forma fechada',
                  'O integrador não encontrou uma primitiva fechada. A integral formal permanece indicada.')
            else:
                output = _result(expr, answer, 'Simbólico',
                  'Primitiva simbólica. Some a constante de integração C.')
                # Uma identidade algébrica nula permite indicar verificação, sem inferir além dela.
                try:
                    output['verified'] = bool(sp.simplify(sp.diff(answer, _X) - expr) == 0)
                except Exception:
                    output['verified'] = False
            output['constant'] = not answer.has(sp.Integral)
            return json.dumps(output, ensure_ascii=False)

        a, b = _parse(low, 'limite inferior'), _parse(high, 'limite superior')
        if a.free_symbols or b.free_symbols:
            raise ValueError('Os limites devem ser números ou ±inf.')
        if (a not in (sp.oo, -sp.oo) and a.is_real is False) or (b not in (sp.oo, -sp.oo) and b.is_real is False):
            raise ValueError('Use limites reais.')
        integral = sp.Integral(expr, (_X, a, b))
        answer = sp.integrate(expr, (_X, a, b))
        if answer in (sp.nan, sp.zoo, sp.oo, -sp.oo):
            return json.dumps(_result(integral, answer, 'Divergente / indefinida',
                'Esta integral não possui valor real finito no sentido impróprio usual; não interpretamos infinito ou nan como uma resposta finita.'), ensure_ascii=False)
        if answer.has(sp.Integral):
            output = _result(integral, answer, 'Integral formal',
                  'Sem avaliação exata conclusiva. Evitamos apresentar uma aproximação como prova de convergência.')
            # Aproximação apenas em intervalo compacto onde a continuidade foi confirmada.
            if a.is_finite and b.is_finite and a.is_number and b.is_number:
                try:
                    interval = sp.Interval(sp.Min(a, b), sp.Max(a, b))
                    domain = continuous_domain(expr, _X, interval)
                    if domain == interval:
                        numeric = sp.N(integral, 15)
                        if numeric.is_real and numeric.is_finite:
                            output['approx'] = str(numeric)
                            output['note'] = 'Estimativa numérica em intervalo de continuidade; não é certificação rigorosa do erro.'
                except Exception:
                    pass
            return json.dumps(output, ensure_ascii=False)
        return json.dumps(_result(integral, answer, 'Exato / simbólico',
            'Resultado simbólico quando resolúvel; divergências e singularidades devem ser interpretadas no domínio real.'), ensure_ascii=False)

    if action == 'limit':
        target = _parse(req.get('target', '0'), 'ponto do limite')
        if target.free_symbols or (target not in (sp.oo, -sp.oo) and target.is_real is False):
            raise ValueError('O ponto do limite precisa ser real ou ±inf.')
        direction = str(req.get('direction', 'both'))
        if direction not in ('both', 'left', 'right'):
            raise ValueError('Direção inválida.')
        if target in (sp.oo, -sp.oo):
            direction = 'right'
        if direction == 'both':
            left = sp.limit(expr, _X, target, dir='-')
            right = sp.limit(expr, _X, target, dir='+')
            same = (left == right)
            if not same:
                try:
                    same = sp.simplify(left - right) == 0
                except Exception:
                    pass
            if not same:
                output = _result(expr, sp.Symbol('DNE'), 'Não existe (bilateral)',
                    'Os limites laterais são diferentes. Não há limite bilateral.',
                    left_tex=_tex(left), right_tex=_tex(right))
                output['approx'] = None
                return json.dumps(output, ensure_ascii=False)
            answer = left
        else:
            answer = sp.limit(expr, _X, target, dir='-' if direction == 'left' else '+')
        if isinstance(answer, sp.AccumBounds) or answer.has(sp.Limit) or answer in (sp.nan, sp.zoo):
            status, note = 'Inconclusivo / não convergente', 'Não foi encontrado um limite pontual finito ou infinito bem definido.'
        else:
            status, note = 'Limite simbólico', 'O resultado usa a análise simbólica de limites, sem estimativa por amostragem.'
        output = _result(expr, answer, status, note)
        output['target_tex'] = _tex(target)
        output['direction'] = direction
        return json.dumps(output, ensure_ascii=False)
    raise ValueError('Operação não reconhecida.')


# ===========================================================================
# Resolução didática: transformações simbólicas verificáveis.
# Uma saída 'steps' contém operações matemáticas efetivas, nunca um roteiro
# genérico chamado indevidamente de demonstração.
# ===========================================================================
def _proof_step(title, tex=None, note=''):
    item = {'title': title, 'note': note}
    if tex:
        item['tex'] = str(tex)
    return item

def _same_symbolic(a, b):
    try:
        return sp.simplify(a - b) == 0
    except Exception:
        return a == b

def _limit_proof(expr, target, direction, output):
    target_tex = _tex(target)
    lim_tex = r'\lim_{x\to ' + target_tex + '}'
    steps = [_proof_step('Escrever o limite', lim_tex + r'\left(' + _tex(expr) + r'\right)',
        'Identificamos a função e o ponto de aproximação.')]
    result_tex = output['result_tex']
    if output.get('status') == 'Não existe (bilateral)':
        steps.append(_proof_step('Comparar os limites laterais',
            r'\lim_{x\to '+target_tex+r'^-}f(x)='+output['left_tex']+
            r',\qquad\lim_{x\to '+target_tex+r'^+}f(x)='+output['right_tex'],
            'Como as tendências laterais são diferentes, o limite bilateral não existe.'))
        return steps

    # Identidade trigonométrica comprovável: 1 - cos(x) = 2 sin^2(x/2).
    if target == 0 and expr.has(sp.cos) and _same_symbolic(expr, (1-sp.cos(_X))/_X**2):
        rewritten = 2*sp.sin(_X/2)**2/_X**2
        if _same_symbolic(expr, rewritten):
            steps.extend([
                _proof_step('Aplicar a identidade do ângulo duplo',
                    r'1-\cos x=2\sin^2\left(\frac{x}{2}\right)',
                    'Usamos cos(x) = 1 − 2sen²(x/2).'),
                _proof_step('Substituir no quociente',
                    lim_tex+r'\frac{2\sin^2(x/2)}{x^2}',
                    'O numerador foi reescrito sem aproximação numérica.'),
                _proof_step('Isolar o limite fundamental',
                    r'\frac12\lim_{x\to0}\left(\frac{\sin(x/2)}{x/2}\right)^2',
                    'Colocamos 1/2 em evidência e fazemos u=x/2.'),
                _proof_step('Usar sen(u)/u → 1',
                    r'\frac12\cdot 1^2=\boxed{\frac12}',
                    'Como u→0, o limite fundamental vale 1.')
            ])
            if result_tex != _tex(sp.Rational(1,2)):
                return steps[:1]  # Sem resultado consistente: não afirmar prova falsa.
            return steps

    # Limites racionais no infinito: divisão pela maior potência.
    if target in (sp.oo,-sp.oo):
        try:
            numerator, denominator = sp.fraction(sp.together(expr))
            pn,pd = sp.Poly(numerator,_X),sp.Poly(denominator,_X)
            if pn.degree() <= 20 and pd.degree() <= 20 and pn.degree() >= 0 and pd.degree() >= 0:
                degree = max(pn.degree(),pd.degree())
                if degree > 0:
                    top = sp.expand(numerator/_X**degree)
                    bottom = sp.expand(denominator/_X**degree)
                    steps.append(_proof_step('Dividir pela maior potência de x',
                        lim_tex+r'\frac{'+_tex(top)+'}{'+_tex(bottom)+'}',
                        'O fator comum x^'+str(degree)+' foi dividido nos dois termos (para x≠0).'))
                    steps.append(_proof_step('Aplicar os limites dos termos de 1/x',
                        r'\lim_{x\to'+target_tex+r'}\frac{'+_tex(top)+'}{'+_tex(bottom)+
                        '}='+result_tex,
                        'Potências 1/x^k desaparecem no infinito; a comparação de graus determina o resultado.'))
                    return steps
        except (sp.PolynomialError,TypeError,ValueError,ZeroDivisionError):
            pass

    # Simplificação algébrica (furos removíveis, fatoração).
    try:
        num,den = sp.fraction(sp.together(expr))
        factored_num=sp.factor(num)
        factored_den=sp.factor(den)
        reduced=sp.cancel(expr)
        reduced_den=sp.fraction(sp.together(reduced))[1]
        true_cancellation = bool(sp.degree(reduced_den,_X) < sp.degree(den,_X))
        if _same_symbolic(reduced,expr) and reduced != expr and true_cancellation:
            if factored_num!=num or factored_den!=den:
                steps.append(_proof_step('Fatorar numerador e denominador',
                    lim_tex+r'\frac{'+_tex(factored_num)+'}{'+_tex(factored_den)+'}',
                    'Procuramos fatores que possam ser cancelados fora do ponto singular.'))
            steps.append(_proof_step('Cancelar fatores comuns',
                lim_tex+r'\left('+_tex(reduced)+r'\right)',
                'O cancelamento vale na vizinhança onde a expressão original está definida.'))
            try:
                direct=sp.limit(reduced,_X,target,dir='+' if direction!='left' else '-')
                if _tex(direct)==result_tex:
                    steps.append(_proof_step('Substituir o ponto após simplificar',
                        r'\boxed{'+result_tex+'}',
                        'O limite da expressão simplificada coincide com o resultado simbólico.'))
                    return steps
            except Exception:
                pass
    except (sp.PolynomialError,TypeError,ValueError,ZeroDivisionError):
        pass

    # L'Hôpital: só mostrar derivações se a condição 0/0 ou inf/inf ocorrer,
    # e se a última expressão for realmente compatível com o resultado.
    if target not in (sp.oo,-sp.oo) and direction in ('both','left','right'):
        try:
            num,den=sp.fraction(sp.together(expr))
            if den.has(_X) and expr.count_ops() <= 70:
                signs=['+','-'] if direction=='both' else ['-' if direction=='left' else '+']
                u,v=num,den
                candidates=[]
                for count in range(1,4):
                    zero_or_inf = True
                    for sign in signs:
                        a,b=sp.limit(u,_X,target,dir=sign),sp.limit(v,_X,target,dir=sign)
                        if not ((a==0 and b==0) or
                            (a in (sp.oo,-sp.oo) and b in (sp.oo,-sp.oo))):
                            zero_or_inf=False
                            break
                    if not zero_or_inf: break
                    du,dv=sp.diff(u,_X),sp.diff(v,_X)
                    if dv == 0: break
                    new_ratio=du/dv
                    candidates.append(_proof_step(
                        'Aplicar L’Hôpital'+(' novamente' if count>1 else ''),
                        lim_tex+r'\frac{'+_tex(du)+'}{'+_tex(dv)+'}',
                        'Numerador e denominador tendem a 0 (ou ambos são infinitos); diferenciamos cada um.'))
                    u,v=du,dv
                    if all(_tex(sp.limit(new_ratio,_X,target,dir=sign))==result_tex
                           for sign in signs):
                        steps.extend(candidates)
                        steps.append(_proof_step('Avaliar o quociente das derivadas',
                            lim_tex+r'\left('+_tex(new_ratio)+r'\right)=\boxed{'+result_tex+'}',
                            'O novo limite coincide com o limite original, sob as hipóteses de L’Hôpital.'))
                        return steps
        except (ValueError,TypeError,NotImplementedError,ZeroDivisionError):
            pass

    # Série de Taylor: apenas caso analítico e conferência simbólica do limite.
    if target.is_finite and expr.count_ops() <= 45:
        try:
            series=sp.series(expr,_X,target,4)
            polynomial=series.removeO() if isinstance(series,sp.Add) and series.has(sp.Order) else None
            if polynomial is not None:
                from_series=sp.limit(polynomial,_X,target,dir='-' if direction=='left' else '+')
                if _tex(from_series)==result_tex:
                    steps.append(_proof_step('Expandir em série de Taylor',
                        _tex(expr)+'='+_tex(series),
                        'A expansão assintótica perto do ponto revela o termo dominante.'))
                    steps.append(_proof_step('Tomar o limite da expansão',
                        r'\lim_{x\to'+target_tex+r'}'+_tex(polynomial)+r'=\boxed{'+result_tex+'}',
                        'Os termos restantes tendem a zero neste caso analítico.'))
                    return steps
        except (ValueError,TypeError,NotImplementedError,ZeroDivisionError):
            pass

    # Substituição direta quando ela realmente coincide com o limite.
    if target.is_finite:
        try:
            direct=sp.simplify(expr.subs(_X,target))
            if direct.is_finite and _tex(direct)==result_tex:
                steps.append(_proof_step('Substituição direta',
                    r'f('+target_tex+r')='+_tex(direct)+r'=\boxed{'+result_tex+'}',
                    'A expressão está definida e é contínua no ponto considerado.'))
                return steps
        except Exception:
            pass

    steps.append(_proof_step('Resultado obtido pelo motor simbólico',
        r'\boxed{'+result_tex+'}',
        'O algoritmo obteve este resultado, mas não foi possível gerar uma cadeia algébrica verificada para esta expressão.'))
    return steps


def _derivative_proof(expr,order,output):
    steps=[_proof_step('Identificar a função e a ordem',
        r'f(x)='+_tex(expr), 'Vamos calcular '+str(order)+' derivada(s) sucessivas.')]
    if expr == _X**2*sp.exp(_X):
        steps.append(_proof_step('Aplicar a regra do produto',
          r'(x^2e^x)^{\prime}=2xe^x+x^2e^x',
          'Tomamos u=x² e v=eˣ: (uv)′ = u′v + uv′.'))
    elif expr.func in (sp.sin,sp.cos,sp.exp) and expr.args[0].has(_X) and expr.args[0]!=_X:
        inner=expr.args[0]
        steps.append(_proof_step('Identificar a função composta',
          r'u='+_tex(inner)+r',\qquad u^\prime='+_tex(sp.diff(inner,_X)),
          'Aplicamos a regra da cadeia: derivada da função externa vezes a derivada do argumento.'))
    previous=expr
    for i in range(1,order+1):
        result=sp.diff(previous,_X)
        steps.append(_proof_step('Derivar '+str(i)+'ª vez',
          r'\frac{d}{dx}\left('+_tex(previous)+r'\right)='+_tex(result),
          'Diferenciamos a expressão do passo anterior.'))
        previous=result
    if _tex(previous)==output['result_tex']:
        steps.append(_proof_step('Confirmar o resultado final',
            r'\boxed{'+_tex(previous)+'}',
            'A derivada de ordem solicitada coincide com a resposta simbólica do SymPy.'))
    if output.get('point_tex') and output.get('point_result_tex'):
        steps.append(_proof_step('Substituir o ponto de avaliação',
            r'f^{('+str(order)+r')}('+output['point_tex']+')='+output['point_result_tex'],
            'Substituição exata na derivada calculada.'))
    return steps[:10]


def _integral_proof(expr,low,high,output):
    definite=bool(low and high)
    integral_tex = r'\int'
    if definite:
        a,b=_parse(low),_parse(high)
        integral_tex+=r'_{'+_tex(a)+r'}^{'+_tex(b)+'}'
    steps=[_proof_step('Escrever a integral',
        integral_tex+_tex(expr)+r'\,dx',
        'Identificamos o integrando e os limites de integração, quando existentes.')]
    if expr.is_polynomial(_X):
        terms=sp.Poly(expr,_X).terms()
        if len(terms)<=8:
            lines=[]
            for (degree,),coef in terms:
                lines.append(_tex(coef*_X**(degree+1)/sp.Integer(degree+1)))
            steps.append(_proof_step('Aplicar a linearidade e a regra da potência',
                r'\int '+_tex(expr)+r'\,dx='+_tex(sp.integrate(expr,_X))+'+C',
                'Integramos termo a termo: ∫xⁿ dx = xⁿ⁺¹/(n+1), para n ≠ −1.'))
    elif expr == _X*sp.exp(_X):
        steps.extend([
            _proof_step('Escolher u e dv',
                r'u=x,\quad dv=e^x\,dx,\quad du=dx,\quad v=e^x',
                'Integração por partes.'),
            _proof_step('Aplicar ∫u dv = uv − ∫v du',
                r'\int xe^x\,dx=xe^x-\int e^x\,dx',
                'Integramos o restante usando a regra exponencial.'),
            _proof_step('Obter a primitiva',
                r'\int xe^x\,dx=(x-1)e^x+C',
                'O resultado fica em forma fechada.')
        ])
    elif expr == _X/(1+_X**2):
        steps.extend([
            _proof_step('Escolher substituição',
                r'u=1+x^2,\quad du=2x\,dx',
                'A derivada do denominador é proporcional ao numerador.'),
            _proof_step('Integrar a fração logarítmica',
                r'\int\frac{x}{1+x^2}\,dx=\frac12\int\frac{du}{u}'+
                r'=\frac12\log(1+x^2)+C',
                'Retornamos a variável x e observamos que 1+x² é positivo.')
        ])
    elif expr == sp.exp(-_X**2):
        steps.append(_proof_step('Reconhecer a função erro',
            r'\int e^{-x^2}\,dx=\frac{\sqrt{\pi}}2\operatorname{erf}(x)+C',
            'A primitiva da gaussiana é expressa pela função especial erf.'))
    elif expr==sp.sin(_X) or expr==sp.cos(_X):
        steps.append(_proof_step('Usar a primitiva trigonométrica',
            r'\int '+_tex(expr)+r'\,dx='+_tex(sp.integrate(expr,_X))+'+C',
            'Aplicamos a derivada inversa do seno ou cosseno.'))

    if definite:
        if output.get('status') in ('Integral formal','Divergente / indefinida'):
            steps.append(_proof_step('Investigar a convergência',
                None,'Este cálculo não apresentou um valor finito conclusivo. É preciso analisar as singularidades e o caráter impróprio da integral.'))
            return steps
        try:
            primitive=sp.integrate(expr,_X)
            if not primitive.has(sp.Integral):
                a,b=_parse(low),_parse(high)
                endpoint=sp.limit(primitive,_X,b,dir='-') if b in (sp.oo,-sp.oo) else primitive.subs(_X,b)
                initial=sp.limit(primitive,_X,a,dir='+') if a in (sp.oo,-sp.oo) else primitive.subs(_X,a)
                if _tex(sp.simplify(endpoint-initial)) == output['result_tex']:
                    steps.append(_proof_step('Aplicar o Teorema Fundamental do Cálculo',
                        r'F(x)='+_tex(primitive)+r',\quad F(b)-F(a)='+
                         _tex(endpoint)+r'-\left('+_tex(initial)+r'\right)',
                        'Avaliamos a primitiva nos extremos, com limites laterais quando necessário.'))
        except Exception:
            pass
        steps.append(_proof_step('Registrar o valor da integral',
            r'\boxed{'+output['result_tex']+'}',
            'Valor simbólico apresentado pelo integrador; confira a convergência em integrais impróprias.'))
    else:
        steps.append(_proof_step('Resultado simbólico do integrador',
            r'\boxed{'+output['result_tex']+'}+C' if output.get('constant') else output['result_tex'],
            'A primitiva é válida em intervalos do domínio onde a derivada recupera o integrando.'))
        if output.get('verified'):
            steps.append(_proof_step('Conferir por derivação',
                r'\frac{d}{dx}\left('+output['result_tex']+r'\right)='+_tex(expr),
                'SymPy simplificou a diferença entre a derivada da primitiva e o integrando a zero.'))
        else:
            steps.append(_proof_step('Ressalva de verificação',
                None,'Não foi reportada prova independente por derivação.'))
    return steps[:10]


# Mantém o cálculo original e adiciona uma resolução simbólica com etapas.
_SCICALC_ORIGINAL_SOLVE=scicalc_pro_solve

def scicalc_pro_solve(request_json):
    raw=_SCICALC_ORIGINAL_SOLVE(request_json)
    output=json.loads(raw)
    if not output.get('ok'):
        return raw
    try:
        req=json.loads(request_json)
        expr=_parse(req.get('expr',''))
        action=req.get('action')
        if action=='limit':
            target=_parse(req.get('target','0'))
            direction=str(req.get('direction','both'))
            output['steps']=_limit_proof(expr,target,direction,output)
        elif action=='derivative':
            output['steps']=_derivative_proof(expr,int(req.get('order',1)),output)
        elif action=='integral':
            output['steps']=_integral_proof(expr,str(req.get('lower','')).strip(),
                         str(req.get('upper','')).strip(),output)
    except Exception:
        # Um erro ao montar a exposição nunca modifica o resultado matemático.
        output['steps']=[_proof_step('Resolução simbólica',
            output.get('result_tex'), 'O motor obteve a resposta, mas não foi possível construir uma cadeia de passos verificados.')]
    return json.dumps(output,ensure_ascii=False)
