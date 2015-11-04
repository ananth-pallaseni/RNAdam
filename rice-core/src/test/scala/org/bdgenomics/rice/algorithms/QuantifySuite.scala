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

  // Create a fake contig fragment from the sequece specified. 
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

  /**
   * Computes the correct results for a run of Index upon a set of fragments described by the 
   * sequences and transcripts given.
   *
   * @param sequences A list of Strings which describe the base sequences of the contig fragments we are 
   *                  emulating
   * @param transcriptNames A list of the names of the transcripts corresponding to the sequences in the 
   *                    previous parameter.
   * @return Returns A manually computed index result
   */
  def expectedResults(sequences: Array[String], transcriptNames: Array[String], kmerLength: Int = 16): Map[(Long, Boolean, String), Map[String, Long]] = {
    val combined = sequences.zip(transcriptNames) // [ (sequence, transcriptName)]
    val intmers = combined.map(tup => (IntMer.fromSequence(tup._1), tup._2)) // [ ( [Intmers] , transcriptName ) ]
    val kmers = intmers.flatMap(tup => { tup._1.map(imer => ((imer.longHash, imer.isOriginal, {if(imer.isOriginal) imer.toCanonicalString else imer.toAntiCanonicalString}), tup._2) )} ) //[ ((hash, orig), transcriptName) ]
    val partialMap = kmers.groupBy(_._1) // Map[ (hash, orig) -> [ ((hash, orig), transcripts)] ]
    val idx = partialMap.mapValues(v => v.groupBy(w => w._2).mapValues(_.size.toLong)) // Map[ (hash, orig) -> Map[ transcriptName -> count ] ]
    idx
  }

  // Compute actual Index and the expected results for Index based on the input sequences
  def createIndex(sequences: Array[String]): (Map[(Long, Boolean), Map[String, Long]], Map[(Long, Boolean, String), Map[String, Long]]) = {
    // Name all the sequences in order
    val names = {for (i <- 0 to sequences.size) yield "seq" + i.toString}.toArray

    // Create a one contig fragment per sequence
    val frags = sc.parallelize( sequences.zip(names).map(c => createContigFragment(c._1, c._2)) )

    // Create an arbitrary set of transcripts (only required as an ancillary argument to Index)
    val transcripts = sc.parallelize( Seq(Transcript("one", Seq("one"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]()),
      Transcript("two", Seq("two"), "gene1", true, Iterable[Exon](), Iterable[CDS](), Iterable[UTR]())) )

    // Compute Index
    val (imap, tmap) = Index(frags, transcripts)

    // Compute expected results
    val expected = expectedResults(sequences, names)

    (imap, expected)
  }

  // Compare the results of Index with expected results
  def compareResults(recieved: Map[(Long, Boolean), Map[String, Long]], expected: Map[(Long, Boolean, String), Map[String, Long]]) : Boolean = {
    // Make a set of the kmers in Expected that matches the format of the kmers in Recieved:
    val formattedExpected = expected.keySet.map(k => (k._1, k._2))

    // Check Sizes
    val recievedSize = recieved.size 
    val expectedSize = expected.size
    val sizesMatch = recievedSize == expectedSize
    // Edge case if either Recieved or Expected is empty
    val sizeMsg = "Recieved Index of Size: " + recievedSize.toString + ", but expected size of " + expectedSize.toString
    if (recievedSize == 0 || expectedSize == 0) {
      println("Index size was zero")
      println(sizeMsg)
      return false
    }

    // Check if all kmers in expected are present in recieved 
    val kmersPresent = recieved.keySet == formattedExpected

    // Check if there are kmers missing
    val missingKmers = expected.keySet.filter(k => {
      val key = (k._1, k._2)
      val inRecieved = recieved.keySet.contains(key)
      !inRecieved
    })
    val missingMsg = if (missingKmers.size > 0) "Kmers missing from Index:\n" + missingKmers.map(k => k.toString + "\n").reduce(_+_) else "No missing kmers"

    // Check if there were kmers added
    val addedKmers = recieved.keySet.filter(k => {
      val inExpected = formattedExpected.contains(k)
      !inExpected
    })
    val addedMsg = if (addedKmers.size > 0) "Kmers added to Index (that shouldn't be present):\n" + addedKmers.map(k => k.toString + "\n").reduce(_+_) else "No added kmers"

    val kmersCorrect = sizesMatch && kmersPresent
    val kmerMsg = sizeMsg + "\n" + missingMsg + "\n" + addedMsg + "\n"
    /*val correct = sizesMatch && kmersPresent
    if (!correct) {
      println("Index was incorrect")
      println(sizeMsg)
      println(missingMsg)
      println(addedMsg)
      println("All Kmers in Recieved Index:")
      recieved.foreach(println(_))
      println("\n")
      println("All expected Kmers")
      expected.foreach(println(_))
    }

    correct*/

    // Check if counts on the transcripts are correct:
    val transcriptsCheck = expected.keySet.map(e => {
      // For each expected kmer, check if it exists in recieved
      val key = (e._1, e._2)
      val keyInRecieved = recieved.contains(key)
      if (keyInRecieved) {
        // If the kmer does exist, then check if the expected transcripts match the recieved ones
        val recMap = recieved.get(key)
        val expMap = expected.get(e)
        val transcriptsCorrect = expMap.keySet.map(t => {
          val transcriptInRecieved = recMap.contains(t)
          if (transcriptInRecieved) {
            // If the transcript exists, then check if the count is accurate
            val recCount = recMap.get(t)
            val expCount = expMap.get(t)
            if (recCount == expCount) {
              (true, "")
            }
            else {
              val transcriptsMsg = "Count for transcript=" + t.toString + " for kmer=" + e.toString + " was " recCount + ". It should be " + expCount.toString
              (false, transcriptsMsg)
            }
          }
          else {
            val transcriptsMsg = t.toString + " not part of Map for " + e.toString + ". It should be."
            (false, transcriptsMsg)
          }
          })
        transcriptsCorrect.map(t => (t._1, t._2 + "\n")).reduce((t1, t2) => (t1._1 && t2._1, t1._2 + t2._2))
      }
      else {
        val transcriptsMsg = e.toString + " not in Index. It should be."
        (false, transcriptsMsg)
      }
      })
    val (transcriptsCountsCorrect, transcriptsMsg) = transcriptsCheck.map(t => (t._1, t._2 + "\n")).reduce((t1, t2) => (t1._1 && t2._1, t1._2 + t2._2))

    val correct = kmersCorrect && transcriptCountsCorrect
    val msg = kmerMsg + "\n" + transcriptsMsg
    if (correct) {
      true
    }
    else {
      println(msg)
      false
    }

    
  }

  // For a given set of input sequences, test if Index produces the correct output
  def testOfIndex(sequences: Array[String]): Boolean = {
    val (recieved, expected) = createIndex(sequences)
    val correct = compareResults(recieved, expected)
    correct
  }

  sparkTest("Simple Test of Index") { 
    /*val testSeq = "ACACTGTGGGTACACTACGAGA"
    val (imap, tmap) = createTestIndex(testSeq)

    // Test kmer mapping
    val imers = testSeq.sliding(16).map(s => IntMer(s))

    assert(imap.size == 7) // 7 kmers of length 16
    
    assert(imers.forall(i => imap((i.longHash, i.isOriginal))("ctg") == 1))


    // Test transcript mapping
    assert(tmap("one").id == "one")
    assert(tmap("two").id == "two")*/
  }

  sparkTest("Less Simple Test of Index") {
    val seq1 = "AAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAA"
    val seq2 = "AAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTT"
    val testPassed = testOfIndex(Array(seq1, seq2))
    assert(testPassed)
  }

  sparkTest("Index - Test to make sure kmers with Canonicality=False make it through") {
    val seq1 = "TTTTTTTTTTTTTTTT"
    val testPassed = testOfIndex(Array(seq1))
    assert(testPassed)
  }

  /*sparkTest("Simple Test of Mapper") {
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
  }*/
  
}
