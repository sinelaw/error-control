%if False
\begin{code}
{-- LANGUAGE NoMonomorphismRestriction --}

module ErrorControl where

import qualified Data.Packed.Vector as Vector
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

\newtheorem{example}{Example}
\newtheorem{definition}{Definition}

\begin{document}

\title{Error Control Coding using Haskell}
\author{Noam~Lewis [lenoam at gmail.com]}

\maketitle

\begin{abstract}
  Error control codes (ECCs) are used to protect data from errors. The theory of ECCs involves such mathematics as finite fields and linear algebra. This report describes an attempt to implement a few ECCs in Haskell. The code takes advantage of the HMatrix package and others. 
\end{abstract}

\section{Introduction}

A Galois field of size $q$ is a finite field with $q$ distinct elements and is denoted $GF(q)$. Our scope will be limited to the case where $q$ is prime. In this case, $GF(q)$ is equivalent to the integers with arithmetic modulo $q$. In most of the cases we'll consider, $q=2$.

\section{Linear Block Codes}

\begin{definition}[Code word]
  A code word of size $n$ is an $n$-tuple of symbols, $(c_0, c_1, \ldots, c_{n-1})$. If the individual elements belong to a Galois field of size $q$, the code word is called $q$-ary.
\end{definition}

The input data is a stream of symbols from $GF(q)$. If the data is sliced into blocks of size $k$, the set of all possible data blocks, that is - all $k$-tuples $(m_0, m_1, \ldots, m_{k-1})$, forms a vector space over $GF(q)$. In that case, there are $q^k$ possible $k$-symbol data vectors. 

\begin{definition}[Block ECC]
  A block error control code $\mathbf{C}$ is a function from input data vectors, to $n$-symbol code words of output data. The number of distinct outputs is denoted $M$. The code words of an ECC are its possible outputs.
\end{definition}

If $M=q^k$ then the input vectors can be of equal length, $k$. This case, of equal input block length (where indeed $M=q^k$) is the only case we shall consider here. The inputs are thus vectors of size $k$, and the outputs are coded vectors of size $n$. Since there are $q^n$ coded words, if $q^n > M$ (which in our case, translates to $n > k$), then the code is redundant. 

\begin{definition}[Code redundancy]
  The redundancy of a code of size $M=q^k$ and output words of length $n$, is defined as $r = n-k$.
\end{definition}


\subsection{Additive noise channel}
The whole point of ECCs is that they can be used to detect and possibly correct errors in the data. The model considered here is the \emph{additive noisy channel model}. Recall that the code words we consider are formed of symbols from $GF(q)$. The additive model then assumes that the original encoding of the data, the \emph{transmitted codewords}, pass through a channel that adds error symbols (via modulo-$q$ addition) to arbitrary places in the codeword. Thus, if the transmitted word is $\mathbf{c}=(c_0, c_1, \ldots, c_{n-1})$, and the error word happens to be $\mathbf{e}=(e_0,e_1,\ldots,e_{n-1})$, then the received word is $\mathbf{r}=\mathbf{c}+\mathbf{e}=(c_0+e_0, c_1+e_1, \ldots, c_{n-1}+e_{n-1})$. Hopefully, most of the $e_i$ symbols are zero. The vector $\mathbf{e}$ is called an \emph{error pattern}. We can thus write: $\mathbf{e} = \mathbf{r} - \mathbf{c}$, or in Haskell:
\begin{code}
  errorPat r c = r - c
\end{code}
Which can be reduced to:
\begin{code}
  errorPat' = (-)
\end{code}


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
weightV = Vector.foldVector (\x -> if x == 0 then id else succ) 0
\end{code}
Here and henceforth, ``Vector'' is the qualified name for HMatrix's ``Data.Packed.Vector''. In many cases we want a measure for the ``distance'' between two words. The most common measure used is the Hamming distance.

\begin{definition}[Hamming distance]
  The Hamming distance between two vectors is the number of elements which differ between them.
\end{definition}

In Haskell, we can implement the Hamming distance function directly, using zipVector (the equivalent to zipWith for Vectors):
\begin{code}
if' a b p = if p then a else b
hammingDist = Vector.foldVector (if' succ id) . Vector.zipVector (==) 
\end{code}
In functional programming we have the luxury of equational reasoning. In the last definition we can reason that elements in vectors are equal if the difference between them is zero. (In other words, $x=y \Leftrightarrow x-y=0$). Thus,
\begin{code}
hammingDist2 = Vector.foldVector (if' succ id) . Vector.zipVector (\x y -> x-y == 0) 
\end{code}
Moving the condition $=0$ into the fold yields:
\begin{code}
hammingDist3 = Vector.foldVector (\e -> if e == 0 then id else succ) . Vector.zipVector (\x y -> x-y) 
\end{code}
Now we notice that the fold matches exactly the above definition for weightV.
\begin{code}
hammingDist4 = weightV . Vector.zipVector (\x y -> x-y) 
\end{code}
Finally, the zip is simply pointwise subtraction, so we have:
\begin{code}
hammingDist5 = weightV . (-)
\end{code}
Recall that the function errorPat' was defined as (-). Therefore,
\begin{code}
hammingDist' = weightV . errorPat'
\end{code}
This last implementation shows that the Hamming distance between two vectors is equivalent to the weight of the error pattern between them. This equivalnce is not hard to see from the start, so the reasoning sequence above may seem redundant. It was presented as a nice example of equational reasoning in Haskell.

\begin{definition}[Minimum Distance of a Block Code]
  The minimum distance of a block code $\mathbf{C}$ is the minimum Hamming distance between every two distinct code words in $\mathbf{C}$.
\end{definition}
If a block code has a minimum distance of $d_{min}$, a code word can only be confused with another if at least $d_{min}$ errors are introduced. Otherwise, there is no way for it to resemble another code word. Thus, the block code can detect all errors that have resulted from less than $d_{min}$ symbol modifications. Another way to put this, is that all error patterns of weight less than $d_{min}$ are detectable. Note that this condition is a global threshold for arbitrary block codes to detect \emph{all} error patterns of a given weight. Many good codes can detect some error patterns with weight $d_{min}$ or more.

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



\end{document}