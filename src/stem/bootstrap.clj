(ns stem.bootstrap
  (:use [stem.constants])
  (:require [stem.util :as u]
            [stem.gene-tree :as gt]
            [clojure.java.shell :as shell]
            [clojure.pprint :only cl-format :as pp]))


(defrecord  DNASeq [name dna-seq])
(defrecord PhylipFile [nseqs nsites dna-seqs])


(defn sample-dna-matrix
  "Given a vector of vectors of dna sequences, returns a bootstrap sample
  Instead of returning another matrix in vector-vector form, it just returns
  a lazy seq of the rows concatenated, since all that's left is to write
  to a file."
  ([m] (sample-dna-matrix m default-seed))
  ([m seed]
     (let [nsites (count (first m))
           rg (u/rand-generator seed)
           rand-col-idxs (repeatedly nsites #(rg :int nsites))] ; samples cols
       (for [x (range (count m)) y rand-col-idxs]
         (get-in m [x y])))))

;; functions for parsing PHYLIP files
(defn print-matrix [m]
  (doseq [v m] (println v)))

(defn str->chars [s]
  (into (vector-of :char) s))

(defn chars->str [v]
  (apply str v))

(defn dna-seqs->dna-matrix
  [dna-seqs]
  (into [] (map :dna-seq dna-seqs)))

(defn phylip-file-parts->dna-seqs
  [s]
  (map
   (fn [[name dna-str]] (DNASeq. name (str->chars dna-str)))
   (partition 2 s)))

(defn file->phylip
  [f]
  (let [file-parts (->> f slurp (.split #"\s+"))
        [nseqs nsites & data] file-parts
        dna-seqs (phylip-file-parts->dna-seqs data)]
    (PhylipFile. (Integer/parseInt nseqs) (Integer/parseInt nsites) dna-seqs)))

(defn write-ssa-infile
  [nseqs nsites names sample]
  (let [header (str nseqs " "  nsites " 1 0 0 0 0")
        ;pads name if not 10 chars - ssa requirement
        pad-name #(pp/cl-format nil "~10A" %)
        dna-seqs (map #(str (pad-name %1) (chars->str %2)) names (partition nsites sample))]
    (u/write-lines
     "infile"
     (cons header dna-seqs))))

(def ssa-settings-file-default
     "include_gaps: 0
burnin: 100
pburnin: 0.05
pbound_un: 0.05
cbound: 10
print_data_ind: 0
print_data_phylip: 0
cbeta: 0.25
abeta: 0.5
percent_mu: 0.25
num_saved_trees: 1
bound_total_iter: 500000
nj_lik: -891.85
est_Kprime: 1
lin_rates: 0
P_accept: 0.90
delta_L_accept: 10.0
length_est: 0
num_boot: 100
seedj: 12345
seedk: 56789")

(defn write-ssa-settings-file
  [model]
  (u/write-to-file
   "settings-ssa"
   (str "model: " model "\n" ssa-settings-file-default)))

(defn parse-ssa-treefile
  [file theta]
  (let [[rate-str _ newick-str] (.split (slurp file) "\n\n")
        rate (first (.split rate-str " "))]
    (gt/parse-gene-tree-str (str rate newick-str) theta)))

(defn run-ssa [exe]
  (shell/sh (str "resources/" exe)))

(defn clean-up []
  (u/delete-files "settings-ssa" "infile" "besttree" "treefile" "paupfile.nex" "results"))

(defn phylip->genetree
  ([p model theta]
     (let [sample (->> (dna-seqs->dna-matrix (:dna-seqs p))
                       (sample-dna-matrix))]
       (phylip->genetree p model theta sample)))
  ([p model theta sample]
     (let [specs (map :name (:dna-seqs p))]
       (try
         (do
           (write-ssa-infile (:nseqs p) (:nsites p) specs sample)
           (write-ssa-settings-file model)
           (run-ssa u/ssa-for-os)
           (parse-ssa-treefile "treefile" theta))
         (catch Exception e
           (clean-up)
           (u/abort "An error occured generating the bootstrap samples...\n" e))
         (finally
          (clean-up))))))

(defn phylip-file->genetree
  [f model theta]
  (phylip->genetree (file->phylip f) model theta))
