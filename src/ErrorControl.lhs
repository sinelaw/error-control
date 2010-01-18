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
The whole point of ECCs is that they can be used to detect and possibly correct errors in the data. The model considered here is the \emph{additive noisy channel model}. Recall that the code words we consider are formed of symbols from $GF(q)$. The additive model then assumes that the original encoding of the data, the \emph{transmitted codewords}, pass through a channel that adds error symbols (via modulo-$q$ addition) to arbitrary places in the codeword. Thus, if the transmitted word is $\mathbf{c}=(c_0, c_1, \ldots, c_{n-1})$, and the error word happens to be $\mathbf{e}=(e_0,e_1,\ldots,e_{n-1})$, then the received word is $\mathbf{r}=\mathbf{c}+\mathbf{e}=(c_0+e_0, c_1+e_1, \ldots, c_{n-1}+e_{n-1})$. Hopefully, most of the $e_i$ symbols are zero. The vector $\mathbf{e}$ is called an \emph{error pattern}.

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

Or, using the weight function on the difference between Vectors:
\begin{code}
hammingDist' = weightV . (-)
\end{code}

\begin{definition}[Minimum Distance of a Block Code]
  The minimum distance of a block code $\mathbf{C}$ is the minimum Hamming distance between every two distinct code words in $\mathbf{C}$.
\end{definition}
If a block code has a minimum distance of $d_{min}$, a code word can only be confused with another if at least $d_{min}$ errors are introduced. Otherwise, there is no way for it to resemble another code word. Thus, the block code can detect all errors that have resulted from less than $d_{min}$ symbol modifications. Another way to put this, is that all error patterns of weight less than $d_{min}$ are detectable. Note that this condition is a global threshold for arbitrary block codes to detect \emph{all} error patterns of a given weight. Many good codes can detect some error patterns with weight $d_{min}$ or more.

\subsection{Minimizing decoder errors}
If the decoder receives an erroneous word, it must somehow pick a valid code word by some criteria. The ``picking'' problem is sometimes known as a classification problem. A natural choice for classification criteria is to try to minimize the probability of error. Denote the following:
\begin{itemize}
\item $p_t(\mathbf{c})$ = probability of $\mathbf{c}$ to be the original (sent) code word;
\item $p_r(\mathbf{r})$ = probability of $\mathbf{r}$ to be the received code word;
\item $\mathbf{c_j}$ = the actual code word sent, $0 \le j<M$ where $M$ is the number of possible distinct code words (the size of the block code);
\item $\mathbf{\bar{c}}$ = the code word picked after the error correction attempt.
\end{itemize}
What we want is to find the code word $c_i$ that minimizes error:
\begin{equation*}
  P(\mbox{error}) = P(\mathbf{\bar{c}} \neq \mathbf{c_j}) = P\left(\bigcup_{i=0,i\neq j}^M{ ( \mathbf{\bar{c}} = \mathbf{c_i} ) }\right) = \sum_{i=0,i \neq j}^M{P(\mathbf{\bar{c}} = \mathbf{c_i})}.
\end{equation*}


\end{document}