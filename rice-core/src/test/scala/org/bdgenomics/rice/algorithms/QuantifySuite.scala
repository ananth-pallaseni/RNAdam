/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.bdgenomics.rice.algorithms

import org.apache.spark.rdd.RDD
import org.bdgenomics.formats.avro.{ Contig, NucleotideContigFragment, AlignmentRecord }
import org.bdgenomics.adam.models.{ Exon, ReferenceRegion, Transcript, CDS, UTR }
import org.bdgenomics.rice.models.{ KmerIndex, IndexMap }
import org.bdgenomics.rice.algorithms.alignment.{ AlignmentModel, SimpleAlignmentModel }
import org.bdgenomics.rice.utils.riceFunSuite
import net.fnothaft.ananas.models._
import net.fnothaft.ananas.avro.{ Kmer, Backing }
import net.fnothaft.ananas.debruijn.ColoredDeBruijnGraph

object TestAlignmentModel extends AlignmentModel {

  def processRead(iter: Iterator[CanonicalKmer],
                  kmerIndex: KmerIndex): Map[String, Double] = {

    iter.toList
      .flatMap(c => kmerIndex.getTranscripts(c))
      .groupBy(_._1) // Map[ TranscriptID -> List[(TranscriptId, Occurrences)] ] 
      .map(v => v._2.reduce((a, b) => (a._1, a._2 + b._2)))
      .map(v => (v._1, v._2.toDouble))
  }
}

class QuantifySuite extends riceFunSuite {

  // Blatant copy of function from ananas/models/ContigFragment to get around protected status
  def buildFromNCF(fragment: NucleotideContigFragment): ContigFragment = {
    // is this the last fragment in a contig?
    val isLast = fragment.getFragmentNumber == (fragment.getNumberOfFragmentsInContig - 1)

    val sequence = IntMer.fromSequence(fragment.getFragmentSequence)
      .map(_.asInstanceOf[CanonicalKmer])

    new ContigFragment(fragment.getContig.getContigName,
      sequence,
      isLast,
      Option(fragment.getFragmentStartPosition).fold(0)(_.toInt))
  }

  def createContigFragment(sequence: String, name: String): ContigFragment = {
    val ncf = NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName(name)
        .build())
      .setFragmentNumber(0)
      .setNumberOfFragmentsInContig(1)
      .setFragmentSequence(sequence)
      .build()

    buildFromNCF(ncf)
  }

  def createTestIndex(sequence: String = "ACACTGTGGGTACACTACGAGA"): (Map[(Long, Boolean), Map[String, Long]], Map[String, Transcript]) = {
    val ncf = NucleotideContigFragment.newBuilder()
      .setContig(Contig.newBuilder()
        .setContigName("ctg")
        .build())
      .setFragmentNumber(0)
      .setNumberOfFragmentsInContig(1)
      .setFragmentSequence(sequence)
      .build()

    val frag = sc.parallelize(Seq(buildFromNCF(ncf)))

    val tx = Seq(Transcript("one", Seq("one"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()),
      Transcript("two", Seq("two"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()))
    val transcripts = sc.parallelize(tx)

    Index(frag, transcripts)
  }

  sparkTest("Simple Test of Index") {
    val testSeq = "ACACTGTGGGTACACTACGAGA"
    val (imap, tmap) = createTestIndex(testSeq)

    // Test kmer mapping
    val imers = testSeq.sliding(16).map(s => IntMer(s))

    assert(imap.size == 7) // 7 kmers of length 16
    imap.foreach(println)
    assert(imers.forall(i => imap((i.longHash, i.isOriginal))("ctg") == 1))

    // Test transcript mapping
    assert(tmap("one").id == "one")
    assert(tmap("two").id == "two")
  }

  sparkTest("Less Simple Test of Index") {
    // Two sequences with repeats of kmer AAAAAAAAAAAAAAAA
    val seq1 = "AAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAA"
    val seq2 = "AAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTT"
    val name1 = "seq1"
    val name2 = "seq2"
    val frags = sc.parallelize(Seq(createContigFragment(seq1, name1), createContigFragment(seq2, name2)))

    val tx = Seq(Transcript("one", Seq("one"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()),
      Transcript("two", Seq("two"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()))
    val transcripts = sc.parallelize(tx)

    val repeat = "AAAAAAAAAAAAAAAA"
    val repeatHash = IntMer(repeat).longHash
    val repeatOriginality = IntMer(repeat).isOriginal
    println("%s -> %#016x <- %s".format(repeat, repeatHash, repeatOriginality))
    val seq1Hashes = IntMer.fromSequence(seq1).map(i => i.longHash)
    val seq2Hashes = IntMer.fromSequence(seq2).map(i => i.longHash)

    val (imap, tmap) = Index(frags, transcripts)

    // With a kmer length of 16, we should have 33 + 18 kmers of which 4 are repeats of AAAAAAAAAAAAAAAA
    imap.foreach(kv => println("%#016x, %s -> %s".format(kv._1._1, kv._1._2, kv._2)))
    assert(imap.size == 33 + 18 - 4 + 1)
    imap.foreach(v => {
      val ((kHash, kStrand), kMap) = v
      println("on " + v)
      // Assert that the kmer exists in the sequences we have presented
      assert(seq1Hashes.contains(kHash) || seq2Hashes.contains(kHash))

      // Check if the kmer is the repeat
      if (kHash == repeatHash && kStrand == repeatOriginality) {
        // Should have two values (one for each sequence the repeat was found in)
        assert(kMap.size == 2)

        // The two values should correspond to the two sequence names
        val names = kMap.map(m => m._1).toList
        assert(names.contains(name1) && names.contains(name2))

        // Each sequence should have two occurrences each
        assert(kMap(name1) == 2 && kMap(name2) == 2)
      } else {
        // If the kmer is not the repeat, then it should only have one value in map
        assert(kMap.size == 1)

        // The value should correspond to one of the two sequence names
        val names = kMap.toList(0)
        assert(names._1 == name1 || names._1 == name2)

        // The sequence should have one occurrence
        assert(names._2 == 1)
      }
    })
  }

  /*  sparkTest("Simple Test of Mapper") {
    val testSeq = "ACACTGTGGGTACACTACGAGA"
    val ar = Array({
      AlignmentRecord.newBuilder()
        .setSequence(testSeq)
        .build()
    })

    val reads = sc.parallelize(ar)

    val (imap, tmap) = createTestIndex(testSeq)

    val kmerIndex = IndexMap(16, imap)

    val m = Mapper(reads, kmerIndex, TestAlignmentModel).collect()

    // Only one read, so only one item in the map
    assert(m.size == 1)

    // The one item should map to ( "ctg" -> number of occurrences of all kmers in contig ctg = 7, as all 7 kmers are from the read)
    assert(m.forall(v => v._2("ctg") == 7.toDouble))

  }

  sparkTest("Less Simple Test of Mapper") {
    val testSeq1 = "ACACTGTGGGTACACTACGAGA"
    val testSeq2 = "CCAGTGACTGGAAAAA"
    val ar = Array({
      AlignmentRecord.newBuilder()
        .setSequence(testSeq1)
        .build()
    },
      {
        AlignmentRecord.newBuilder()
          .setSequence(testSeq2)
          .build()
      })

    val reads = sc.parallelize(ar)

    val (imap, tmap) = createTestIndex(testSeq1 + testSeq2)

    val kmerIndex = IndexMap(16, imap)

    val m = Mapper(reads, kmerIndex, TestAlignmentModel).collect()

    // Two reads, so two items in the map
    assert(m.size == 2)

    // Assert that readIDs are unique
    assert(m(0)._1 != m(1)._1)

    // Should have either 7 or 1 kmers per read (all kmers are from the same transcript "ctg")
    assert(m.forall(v => { v._2("ctg") == 7.toDouble || v._2("ctg") == 1 }))

    // Should record 8 kmers in total
    assert(m(0)._2("ctg") + m(1)._2("ctg") == 8)

  }

  sparkTest("Testing SimpleAlignmentModel") {
    val testSeq = "ACACTGTGGGTACACTACGAGA"
    val ar = Array({
      AlignmentRecord.newBuilder()
        .setSequence(testSeq)
        .build()
    })

    val reads = sc.parallelize(ar)

    val (imap, tmap) = createTestIndex(testSeq)

    val kmerIndex = IndexMap(16, imap)

    val m = Mapper(reads, kmerIndex, SimpleAlignmentModel).collect()

    // All 7 kmers in this read came from transcript "Ctg", readLength = 22, kmerLength = 16
    // so likelihood = 7 / (22 - 16 + 1) = 1
    assert(m(0)._2("ctg") == 1D)
  }
  */
}
