%if False
\begin{code}
{-# LANGUAGE QuasiQuotes, NoMonomorphismRestriction #-}

module ErrorControl where

import qualified Data.Packed.Static as Static
\end{code}
%endif

\documentclass[onecolumn,x11names,twoside,a4paper,english]{IEEEtran}
\usepackage[english]{babel}
\usepackage[pdftex]{graphicx}
\usepackage{amssymb}
\usepackage{amsmath}
\usepackage{caption}
\usepackage{float}
\usepackage{tikz}
\usepackage{euler}                                %Nicer numbers
\usepackage{subfigure}
\usepackage[shell,pdf]{dottex}
%include polycode.fmt 

\newtheorem{theorem}{Theorem}
\newtheorem{example}{Example}
\newtheorem{definition}{Definition}

\begin{document}

\title{Error Control Coding using Haskell}
\author{Noam~Lewis [lenoam at gmail.com]}

\maketitle

\begin{abstract}
  Error control codes (ECCs) are used to protect data from errors. The theory of ECCs involves such mathematics as finite fields and linear algebra. This report describes an attempt to implement a few ECCs in Haskell. The code takes advantage of the hmatrix-static package and others. 
\end{abstract}

\section{Introduction}

A Galois field of size $q$ is a finite field with $q$ distinct elements and is denoted $GF(q)$. A theorem states that $q$ must be a power of a prime, $q=p^r$ where $p$ is prime. If $q$ is prime, $GF(q)$ is equivalent to the integers with arithmetic modulo $q$. 

\section{Linear Block Codes}

A block code takes vectors (blocks) of input symbols and encodes each one into an output vector (code word). The input data is a stream of symbols from $GF(q)$. If the data is sliced into blocks of size $k$, the set of all possible data blocks, that is - all $k$-tuples $(m_0, m_1, \ldots, m_{k-1})$, forms a vector space over $GF(q)$. In that case, there are $q^k$ possible $k$-symbol data vectors. 

\begin{definition}[Block ECC]
  A block error control code $\mathbf{C}$ is a function from input data vectors, to $n$-symbol code words of output data. The number of distinct outputs is denoted $M$. The code words of an ECC are its possible outputs.
\end{definition}

\begin{definition}[Code word]
  A code word is a value in the image (codomain) of a block ECC. A code word of size $n$ is an $n$-tuple (a vector) of symbols, $(c_0, c_1, \ldots, c_{n-1})$. If the individual elements belong to a Galois field of size $q$, the code word is called $q$-ary.
\end{definition}

If $M=q^k$ then the input vectors can be of equal length, $k$. This case, of equal input block length (where indeed $M=q^k$) is the only case we shall consider here. The inputs are thus vectors of size $k$, and the outputs are coded vectors of size $n$. Since there are $q^n$ coded words, if $q^n > M$ (which in our case, translates to $n > k$), then the code is redundant. 

\begin{definition}[Code redundancy]
  The redundancy of a code of size $M=q^k$ and output words of length $n$, is defined as $r = n-k$.
\end{definition}

We can translate the above definitions into Haskell:
\begin{code}
  codeSize :: (Integral a) => a -> a -> a
  codeSize q n = q^n

  redundancy :: (Integral a) => a -> a -> a
  redundancy n k = n-k
\end{code}
It's easy to see what the point-less versions of those functions are.

\subsection{Additive noise channel}
The whole point of ECCs is that they can be used to detect and possibly correct errors in the data. The model considered here is the \emph{additive noisy channel model}. Recall that the code words we consider are formed of symbols from $GF(q)$. The additive model then assumes that the original encoding of the data, the \emph{transmitted codewords}, pass through a channel that adds error symbols (via modulo-$q$ addition) to arbitrary places in the codeword. Thus, if the transmitted word is $\mathbf{c}=(c_0, c_1, \ldots, c_{n-1})$, and the error word happens to be $\mathbf{e}=(e_0,e_1,\ldots,e_{n-1})$, then the received word is $\mathbf{r}=\mathbf{c}+\mathbf{e}=(c_0+e_0, c_1+e_1, \ldots, c_{n-1}+e_{n-1})$. Hopefully, most of the $e_i$ symbols are zero. The vector $\mathbf{e}$ is called an \emph{error pattern}. We can thus write: $\mathbf{e} = \mathbf{r} - \mathbf{c}$, or in Haskell:
\begin{code}
  errorPat' :: Vector n a -> Vector n a -> Vector n a
  errorPat' r c = r - c
\end{code}
Which can be reduced to:
\begin{code}
  errorPat :: Vector n a -> Vector n a -> Vector n a
  errorPat = (-)
\end{code}
The statically-sized vector type is provided by hmatrix-static, a thin wrapper around hmatrix. The package uses type-level integers from tfp to allow static checking that vectors match in size.

\subsection{Detecting and reacting to errors}
The decoder recieves encoded words with errors added. Thus, it must first identify whether errors have occured in the received words. This identification is called \emph{error detection}. If the received erroneous word happens to match a different, valid, but unintended code word, then error detection will fail by what we call an \emph{undetectable error}. Since there are $M$ total valid coded words, for each coded word there are $M-1$ other code words that may be mistakenly accepted by the decoder.

Once an error has been detected, the decoder may react in various ways:
\begin{itemize}
\item Throw the erroneous word away, noting that a word is missing. This is useful for real-time lossy communications, such as voice conversations.
\item ``Retry'' a transmission of the erroneous word. In a live communication system, this method is used in Automatic Repeat Request (ARQ) protocols. In other cases, such as reading from physical media, this may amount to attempting another read.
\item Try to correct the errors in the word (\emph{forward error correction}, or FEC).
\end{itemize}

Error correcting systems can make a mistake when attemping to correct an error. This type of mistake is called a \emph{decoder error}. In many cases, some error patterns are more likely to occur than others. The common case is where error patterns have a higher chance of containing mostly non-zero symbols, and the higher the number of non-zero symbols in the pattern, the lower the chance for that error pattern to occur. This assumption is used by FECs to decide how to correct an erroneous code word, and leads to the following definition.

\begin{definition}[Weight]
  The weight of a code word, or error pattern, is the number of non-zero symbols it contains.
\end{definition}

In Haskell, we use Vector-typed values (from the HMatrix package) for input data blocks, code words and error patterns. The weight function can then be defined as:

\begin{code}
weightV :: Vector n a -> Int
weightV = Vector.foldVector (\x -> if x == 0 then id else succ) 0
\end{code}
Here and henceforth, ``Vector'' is the qualified name for HMatrix's ``Data.Packed.Vector''. In many cases we want a measure for the ``distance'' between two words. The most common measure used is the Hamming distance.

\begin{definition}[Hamming distance]
  The Hamming distance between two vectors is the number of elements which differ between them.
\end{definition}

Note that the Hamming distance is commutative. In Haskell, we can implement the Hamming distance function directly, using zipVector (the equivalent to zipWith for Vectors):
\begin{code}
if' :: a -> a -> Bool -> a
if' a b p = if p then a else b

hammDist1 :: Vector n a -> Vector n a -> Int
hammDist1 = Vector.foldVector (if' succ id) . Vector.zipVector (==) 
\end{code}
In functional programming we have the luxury of equational reasoning. In the last definition we can reason that elements in vectors are equal if the difference between them is zero. (In other words, $x=y \Leftrightarrow x-y=0$). Thus,
\begin{code}
hammDist2 = Vector.foldVector (if' succ id) . Vector.zipVector (\x y -> x-y == 0) 
\end{code}
Moving the condition $=0$ into the fold yields:
\begin{code}
hammDist3 = Vector.foldVector (\e -> if e == 0 then id else succ) . Vector.zipVector (\x y -> x-y) 
\end{code}
Now we notice that the fold matches exactly the above definition for weightV.
\begin{code}
hammDist4 = weightV . Vector.zipVector (\x y -> x-y) 
\end{code}
Finally, the zip is simply pointwise subtraction, so we have:
\begin{code}
hammDist5 = weightV . (-)
\end{code}
Recall that the function errorPat' was defined as (-) on vectors. Therefore,
\begin{code}
hammDist :: Vector n a -> Vector n a -> Int
hammDist = weightV . errorPat
\end{code}
This last implementation shows that the Hamming distance between two vectors is equivalent to the weight of the error pattern between them. These conclusions are not hard to see from the start, so the reasoning sequence above may seem redundant. It was presented as a nice example of equational reasoning in Haskell.

\begin{definition}[Minimum Distance of a Block Code]
  The minimum distance of a block code $\mathbf{C}$ is the minimum Hamming distance between every two distinct code words in $\mathbf{C}$.
\end{definition}
If a block code has a minimum distance of $d_{min}$, a code word can only be confused with another if at least $d_{min}$ errors are introduced. Otherwise, there is no way for it to resemble another code word. Thus, the block code can detect all errors that have resulted from less than $d_{min}$ symbol modifications. Another way to put this, is that all error patterns of weight less than $d_{min}$ are detectable. Note that this condition is a global threshold for arbitrary block codes to detect \emph{all} error patterns of a given weight. Many good codes can detect some error patterns with weight $d_{min}$ or more.

Given the set of all code words, we can implement a function that calculates the $d_{min}$. In many we can derive this distance analytically.
\begin{code}
minHammDist' :: [Vector n a] -> Int
minHammDist' (v:vs) = minHammDist vs `min` minDistFrom v vs
    where minDistFrom v = foldr1 min . map (hammDist v)
minHammDist' (v1:v2:[]) = hammDist v1 v2
minHammDist' _      = error "Need two or more vectors to calculate minimum distance"
\end{code}
We can generalize minDist to a function that folds over a commutative binary operation:
\begin{code}
commFoldr :: (b -> b -> b) -> (a -> a -> b) -> [a] -> b
commFoldr mergeF binF (x:xs) = commFoldr mergeF binF xs `mergeF` commFoldr' xs
    where commFoldr' = foldr1 mergeF . map (binF x)
commFoldr _ binF (x1:x2:[]) = binF x1 x2
commFoldr _ _ _ = error "Need two or more values to calculate commutative fold"
\end{code}
Using the generalization simplifies minHammDist' to:
\begin{code}
minHammDist :: [Vector n a] -> Int
minHammDist = commFoldr min hammDist
\end{code}
We can use commFoldr to perform other similar calculations, such as the maximum Euclidean distance between vectors:
\begin{code}
eucDist :: Floating a => Vector n a -> a
eucDist u v = sqrt (w `dot` w)
    where w = u - v

maxEucDist :: [Vector n a] -> a
maxEucDist = commFoldr max eucDist
\end{code}



\subsection{Minimizing decoder errors}
If the decoder receives an erroneous word, it must somehow pick a valid code word by some criteria. The ``picking'' problem is sometimes known as a classification problem. A natural choice for classification criteria is to try to minimize the probability of error. Denote the following:
\begin{itemize}
\item $\mathbf{c_j}$ = the actual code word sent, $0 \le j<M$ where $M$ is the number of possible distinct code words (the size of the block code);
\item $\mathbf{r}$ = the code word received, after additive errors;
\item $p_t(\mathbf{c})$ = probability of $\mathbf{c}$ to be the original (sent) code word;
\item $p_r(\mathbf{r})$ = probability of receiving the code word $\mathbf{r}$;
\end{itemize}
What we want is to find the code word $\mathbf{c_i}$ that minimizes error:
\begin{equation*}
  P(\mbox{picking}\ \mathbf{c_i}\ \mbox{is an error} \mid \mathbf{r}\ \mbox{received}) =  P(\mathbf{c_i} \neq \mathbf{c_j} \mid \mathbf{r}) = 1 - P(\mathbf{c_i} = \mathbf{c_j} \mid \mathbf{r})
\end{equation*}
In other words, to minimize error, we want to \emph{maximize} the probability for correct decoding. The probability we want to maximize is $P(\mathbf{c_i} = \mathbf{c_j} \mid \mathbf{r}) = P(\mathbf{c_i}\ \mbox{was sent} \mid \mathbf{r}\ \mbox{recieved})$, and the criteria is called \emph{maximum a posteriori} decoding. Using Bayes' rule we can transform it into an expression on the \emph{a priori} probability $P(\mathbf{r}\ \mbox{received} \mid  \mathbf{c_i}\ \mbox{sent})$, the probability for receiving $\mathbf{r}$ under the condition that $\mathbf{c_i}$ was sent. This probability can often be calculated from the properties of the channel and the block code used. Using Bayes' rule gives:
\begin{equation*}
  P(\mathbf{c_i} \mid \mathbf{r}) = \frac{p_t(\mathbf{c_i}) P(\mathbf{r} \mid \mathbf{c_i})} {p_r(\mathbf{r})}
\end{equation*}
Recall that we are searching for the code word $\mathbf{c_i}$ that maximizes the above expression. If the probabilities $p_t(\mathbf{c_i})$ are the same for all the code words (i.e. for all $0 \leq i < M$), then we can disregard that as a constant factor that has no effect on the maximum of the right-hand-side expression. In this case, the maximum of the whole expression will coincide with the maximum of $P(\mathbf{r} \mid \mathbf{c_i})$. The constancy condition on $p_t$ is plausible - it means that all code words are equally likely to be transmitted. Picking the code word $\mathbf{c_i}$ that maximizes $P(\mathbf{r} \mid \mathbf{c_i})$ is known as \emph{maximum likelihood} decoding. Under the condition of constant $p_t$, it is equivalent to the aforementioned maximum a posteriori decoding, and it minimizes the probability of error.

We may express $P(\mathbf{r} \mid \mathbf{c_i})$ in error pattern terms. Let $\mathbf{e_i} = \mathbf{r} - \mathbf{c_i}$, be the error pattern corresponding to $\mathbf{c_i}$ having been transmitted, and $\mathbf{r}$ received. Then the probability for receiving $\mathbf{r}$ will be equal to the probability of $\mathbf{e_i}$ having been the error pattern in the transmission. Assuming that lower-weight error patterns are more likely than higher-weight ones, the code word $\mathbf{c_i}$ that \emph{maximizes} $P(\mathbf{r} \mid \mathbf{c_i})$ is the one that corresponds to the \emph{lowest}-weight error pattern. However, the weight of the error pattern $\mathbf{e_i}$ is by definition the number of non-zero symbols in the pattern $\mathbf{r} - \mathbf{c_i}$. This number is exactly the Hamming distance between $\mathbf{r}$ and $\mathbf{c_i}$. To summarize, the minimum error will be achieved if we pick $\mathbf{c_i}$ that is nearest (in Hamming distance) to $\mathbf{r}$.

A consequence of the above statement is that correct decoding depends on receiving a word that is closest to the transmitted code word. The minimum distance between code words is $d_{min}$. Thus, if the received word is more than $d_{min}/2$ away from the sent code word, it may end up closer to a different, unintended code word which will cause the decoder to misclassify. The conclusion is that a block code with minimum distance $d_{min}$ can correct all error patterns of weight less than $d_{min}/2$.

\subsection{Complete and bounded distance decoders}

\begin{definition}[Complete decoder]
  Given a received word $\mathbf{r}$, a complete error correcting decoder picks the closest code word (in Hamming distance) as the result of the decoding.
\end{definition}
For efficiency reasons, we often prefer to use bounded distance decoders, which are defined as:
\begin{definition}[Bounded distance decoder]
  Given a received word $\mathbf{r}$ and a distance bound $t$, a $t$-error correcting bounded distance decoder considers the closest code word (the one with the minimal Hamming distance, $d$) within the sphere of distance $t$. If such a code word is found, it is the result of the decoding. Otherwise, the decoding fails. 
\end{definition}
Note that for the distance bound $t$ to make sense, it must satisfy $t \leq d_{min}/2$. Beyond the distance of $d_{min}/2$ the decoder will always find a code word, and the bound is useless - the decoder would become a complete decoder.

\subsection{Minimal redundancy }
\begin{definition}[Hamming sphere]
  A Hamming sphere of code word $\mathbf{c}$ with radius $t$ contains all the code words that are at a Hamming distance of $t$ or less from $\mathbf{c}$. 
\end{definition}
The volume of a Hamming sphere is the number of code words in the sphere. A known result (supposedly trivial) is that the volume $V_q(n,t)$ of a Hamming sphere, for $q$-ary code words of size $n$ and sphere radius $t$ satifies:
\begin{equation}
  \label{eq:hamm-sphere-vol}
  V_q(n,t) = \sum_{j=0}^t {n \choose j} (q-1)^j
\end{equation}
In the case of $q=2$, the volume is the number of $n$-bit words that differ from the center by $t$ bits or less:
\begin{equation}
  \label{eq:hamm-sphere-vol2}
  V_2(n,t) = \sum_{j=0}^t {n \choose j} = {n+t \choose t}
\end{equation}
In Haskell:
\begin{code}
-- Vectors in a Hamming Sphere of radius t, centered at c
hammSphere :: Integral b => b -> Vector n a -> [Vector n a] -> [Vector n a]
hammSphere t c = filter ((<= t) . (hammDist c))

-- Choice function
choose :: (Integral a) => a -> a -> a
choose n k = product [(n-k+1)..n] `div` product [1..k]

-- Volume of a Hamming sphere
hammSphereVol :: Integral a => a -> a -> a -> a
hammSphereVol 2 n t = ((n+t) `choose` t)
hammSphereVol q n t = sum [(n `choose` j) * (q-1)^j | j <- [0..t]]
\end{code}

Let's say we want a code that can correct $t$ errors. What is the minimal redundancy required? This question is an unsolved problem, but there are known bounds to the answer.
\begin{theorem}[Hamming bound]
  A $q$-ary code with words of length $n$ can correct $t$ errors only if the redundancy $r$ satisfies:
  \begin{equation}
    \label{eq:hamm-bound}
    r \ge log_q{V_q(n,t)}
  \end{equation}
\end{theorem}
This places a lower bound on the required redundancy. There is no point in searching for a code that can correct $t$ but has less redundancy than given by the Hamming bound. There is also an \emph{upper} bound on how much redundancy is required:
\begin{theorem}[Gilbert Bound]
  There exists a $t$-error correcting, $q$-ary code of length $n$ and redundancy $r$ that satisfies:
  \begin{equation}
    \label{eq:gil-bound}
    r \leq log_q{V_q(n,2t)}
  \end{equation}
\end{theorem}
There is no point in using a code that has more redundancy than the Gilbert bound, because there certainly exists a better code. The better code has redundancy equal to the Gilbert bound, and can be constructed by selecting arbitrary code words that are spaced at least $2t$ apart. 

We can summarize the two bounds by the inequality:
\begin{equation*}
    log_q{V_q(n,2t)} \ge r \ge log_q{V_q(n,t)}
\end{equation*}



\subsubsection{Implementing code block generators in Haskell}

The implementation of a code generator (such as the one just mentioned) is problematic, if we want static type checking on the lengths of vectors. For example, if we want to implement a function that generates codes with Gilbert-bound redundancy, we might start with a type signature like:
\begin{code}
maxRedunCode :: Integral a => a -> a -> a -> [Vector n a]
maxRedunCode q m t = -- ??
\end{code}
The function is supposed to generate $m$ code words of length $n$ over $GF(q)$, that can correct $t$ errors. I am unsure how the function can construct the appropriately-sized vectors. One issue is that the hmatrix-static package doesn't seem to offer a way to edit a specific element in the vector.

\subsection{Perfect codes}
\begin{definition}[Perfect code]
  A code is a \emph{perfect code} if it satisfies the equality in the Hamming bound.
\end{definition}
Perfect codes are not the best practical codes. Reed-Solomon codes, for example, are far from perfect in the redundancy sense, but they turn out to be of great practical use. 

There is a theorem that states that for a $q$-ary code to be perfect, it must have size (number of code words) $M=q^k$ for some positive integer $k$. This theorem along with the Hamming bound leads to the following constraint for perfect codes:
\begin{equation}
  \label{eq:hamb-perf}
  \sum_{j=0}^t{{n \choose j}(q-1)^j} = q^{n-k}
\end{equation}
There are five classes of solutions to the above constraint:
\begin{itemize}
\item $n=k,\ t=0$ - These are the trivial codes that have zero redundancy, and cannot detect or correct any errors.
\item $q=2,\ n\ \mbox{odd},\ k=1,\ t=(n-1)/2$ - Odd-length, binary repetition codes. These codes contain two code words of length $n$: all zeros or all ones.
\item $q=2,\ n=2^m-1,\ k=2^m-m-1,\ t=1$ and $m$ is a positive integer - Hamming codes (discussed later) are the only linear codes in this class, but there are also non-linear codes.
\item $q=2,\ n=23,\ k=12,\ t=3$ - Only one code exists in this class, the Golay code $\mathcal{G}_{23}$.
\item $q=3,\ n=11,\ k=6,\ t=2$ - Similarly, only one code exists in this class, the Golay code $\mathcal{G}_{11}$.
\end{itemize}
A theorem due to Tietavainen shows that no perfect codes exist other than those with parameters from the above list. Thus, all other solutions to Equation \ref{eq:hamb-perf} result in parameters that have no corresponding codes.

\begin{theorem}
  Any perfect code must have the same length $n$, symbol field $GF(q)$ and size $M=q^k$ as a Hamming, Golay or repetition code.
\end{theorem}

\end{document}